/*
 * prs_erosion.h — C interfejs biblioteke za hidrauličnu eroziju (PRS-Project).
 *
 * Ovo je JEDINI header koji korisnik DLL-a (npr. opengl-forward-renderer)
 * treba. Namjerno je čisti C:
 *   - nema C++ tipova ni izuzetaka preko granice DLL-a, pa DLL može koristiti
 *     program preveden bilo kojim kompajlerom (MinGW GCC, Clang, MSVC);
 *   - stanje simulacije je neprozirni handle (PrsErosion*), pa se unutrašnjost
 *     (raspored polja, SIMD/OpenMP verzije) može mijenjati bez promjene API-ja.
 *
 * Konvencije mreže (iste kao TerrainHeightfield u rendereru):
 *   - indeks ćelije: y * width + x (red po red); y mreže = osa Z u svijetu
 *   - sve vrijednosti su float, visine u svjetskim jedinicama
 *   - cellSize = razmak između centara susjednih ćelija
 *
 * Model: virtual pipe plitka voda (Mei, Decaudin, Hu 2007) s erozijom,
 * depozicijom, transportom sedimenta i isparavanjem.
 *
 * Stijena (bedrock): uz teren se može zadati visina stijene. Između stijene i
 * površine je sloj tla koji voda može odnijeti; stijenu ne može (ili može
 * sporije, vidi PrsErosionParams::rockErodibility). Gdje je stijena gola,
 * voda teče, a teren se ne mijenja (stalno riječno korito).
 *
 * Funkcije nisu thread-safe za isti handle; različiti handle-ovi su nezavisni.
 */
#ifndef PRS_EROSION_H
#define PRS_EROSION_H

#ifdef __cplusplus
extern "C" {
#endif

/* ── Izvoz simbola ──────────────────────────────────────────────────────── */
#if defined(_WIN32) || defined(__CYGWIN__)
#  if defined(PRS_EROSION_BUILDING_DLL)
#    define PRS_EROSION_API __declspec(dllexport)
#  else
#    define PRS_EROSION_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__)
#  define PRS_EROSION_API __attribute__((visibility("default")))
#else
#  define PRS_EROSION_API
#endif

/* Mijenja se kad se API promijeni na nekompatibilan način. Program može
 * provjeriti da DLL odgovara headeru:
 *   prs_erosion_api_version() == PRS_EROSION_API_VERSION */
#define PRS_EROSION_API_VERSION 2

/* Vrijednost u polju stijene koja znači "u ovoj ćeliji nema stijene"
 * (najmanji konačan float, -FLT_MAX). */
#define PRS_EROSION_NO_BEDROCK (-3.402823466e+38F)

/* Neprozirno stanje jedne simulacije. */
typedef struct PrsErosion PrsErosion;

/* Kako sediment putuje s vodom. */
typedef enum PrsSedimentTransport {
    /* Konzervativno, kroz iste virtualne cijevi kao voda (zadano). */
    PRS_TRANSPORT_FLUX = 0,
    /* Semi-Lagrangian advekcija iz Mei et al.; ne čuva masu (za poređenje). */
    PRS_TRANSPORT_SEMI_LAGRANGIAN = 1
} PrsSedimentTransport;

/* Parametri simulacije. Uvijek prvo popuniti s prs_erosion_default_params(),
 * pa mijenjati samo ono što treba. */
typedef struct PrsErosionParams {
    float dt;                  /* vremenski korak [s] */
    float gravity;             /* g */
    float pipeArea;            /* A, poprečni presjek virtualne cijevi */

    float rainRate;            /* visina vode dodana svakoj ćeliji u sekundi */
    float evaporationRate;     /* Ke: udio vode koji ispari u sekundi */
    float minDepthForVelocity; /* ispod ove dubine brzina = 0 */

    int   erosionEnabled;      /* 0 = samo voda, bez erozije i sedimenta */
    float sedimentCapacity;    /* Kc */
    float dissolvingRate;      /* Ks [1/s] */
    float depositionRate;      /* Kd [1/s] */
    float minTilt;             /* donja granica za sin(alfa) */
    float maxErosionDepth;     /* Kdmax: dubina na kojoj erozija dostiže punu snagu */
    int   sedimentTransport;   /* PrsSedimentTransport */

    /* Tvrdoća stijene u [0, 1]: 0 = stijena se ne erodira (zadano),
     * 0.05 = 20x sporije od tla, 1 = isto kao tlo. Nema uticaja bez stijene. */
    float rockErodibility;
} PrsErosionParams;

/* Postavke scenarija "rijeka" (prs_erosion_setup_river). Prvo popuniti s
 * prs_erosion_default_river_settings(). */
typedef struct PrsRiverSettings {
    float soilDepth;  /* debljina tla izvan korita (u koritu je gola stijena) */
    float sourceRate; /* jačina izvora: visina vode u sekundi po ćeliji izvora */
    float baseHeight; /* visina dna bazena (najniža tačka terena) */
} PrsRiverSettings;

/* Statistika za prikaz u UI-ju. */
typedef struct PrsErosionStats {
    long long steps;         /* broj urađenih koraka od kreiranja / reseta */
    double    simTime;       /* simulirano vrijeme [s] */
    double    lastStepMs;    /* prosječno vrijeme jednog koraka u zadnjem prs_erosion_step [ms] */
    double    totalWater;    /* sum(d) * cellSize^2 */
    double    totalMaterial; /* sum(b + s) * cellSize^2; konstantno uz PRS_TRANSPORT_FLUX */

    long long rockCells;     /* ćelije u kojima je stijena zadana */
    long long bareRockCells; /* od toga gola stijena (tlo tanje od 0.01) */
    double    soilVolume;    /* zapremina tla iznad stijene: sum(b - r) * cellSize^2 */
} PrsErosionStats;

/* ── Verzija i greške ───────────────────────────────────────────────────── */

PRS_EROSION_API int prs_erosion_api_version(void);

/* Opis zadnje greške u ovoj niti (prazan string ako je nije bilo).
 * Funkcije koje mogu pasti vraćaju NULL ili 0 i postavljaju ovu poruku. */
PRS_EROSION_API const char* prs_erosion_last_error(void);

/* ── Kreiranje i uništavanje ────────────────────────────────────────────── */

PRS_EROSION_API void prs_erosion_default_params(PrsErosionParams* out);

/* Nova simulacija na mreži width x height (oba >= 2).
 *   terrain: width*height visina terena; NULL = ravan teren na visini 0
 *   params:  NULL = zadani parametri
 * Voda, sediment, fluks i brzina počinju od nule.
 * Vraća NULL pri grešci (vidi prs_erosion_last_error). */
PRS_EROSION_API PrsErosion* prs_erosion_create(int width, int height, float cellSize,
                                               const float* terrain,
                                               const PrsErosionParams* params);

/* NULL je dozvoljen. */
PRS_EROSION_API void prs_erosion_destroy(PrsErosion* sim);

/* Vraća simulaciju na početak s novim terenom (iste dimenzije): briše vodu,
 * sediment, fluks, brzinu, statistiku i STIJENU (ona zavisi od terena, pa je
 * poslije reseta treba zadati ponovo). Izvori vode i parametri ostaju.
 * terrain NULL = ravan teren na visini 0. Vraća 1 ako uspije. */
PRS_EROSION_API int prs_erosion_reset(PrsErosion* sim, const float* terrain);

/* ── Parametri i izvori vode ────────────────────────────────────────────── */

/* Vraća 1 ako uspije, 0 ako parametri nisu ispravni (npr. dt <= 0). */
PRS_EROSION_API int  prs_erosion_set_params(PrsErosion* sim, const PrsErosionParams* params);
PRS_EROSION_API void prs_erosion_get_params(const PrsErosion* sim, PrsErosionParams* out);

/* Izvor vode: dodaje `rate` visine vode u sekundi svakoj ćeliji unutar kruga
 * poluprečnika `radius` (u ćelijama) oko (x, y). Vraća 1 ako uspije. */
PRS_EROSION_API int  prs_erosion_add_source(PrsErosion* sim, int x, int y, float radius, float rate);
PRS_EROSION_API void prs_erosion_clear_sources(PrsErosion* sim);

/* ── Stijena (bedrock) ──────────────────────────────────────────────────── */

/* Sve funkcije za stijenu koriste TRENUTNI teren i garantuju stijena <= teren.
 * Pozivaju se poslije prs_erosion_create / prs_erosion_reset, prije koraka.
 * Vraćaju 1 ako uspiju. */

/* Mapa visina stijene (width*height vrijednosti). Stijena viša od terena se
 * spušta na teren (gola stijena); NaN ili beskonačna vrijednost znači "ovdje
 * nema stijene". heights == NULL uklanja stijenu. */
PRS_EROSION_API int prs_erosion_set_bedrock(PrsErosion* sim, const float* heights);

/* Debljina tla po ćeliji (width*height vrijednosti): stijena = teren - debljina.
 * 0 ili negativno = gola stijena; NaN ili beskonačno = ovdje nema stijene. */
PRS_EROSION_API int prs_erosion_set_soil_depth(PrsErosion* sim, const float* soilDepth);

/* Isti sloj tla svuda: stijena = teren - soilDepth (prati oblik terena). */
PRS_EROSION_API int prs_erosion_set_uniform_soil(PrsErosion* sim, float soilDepth);

/* ── Scenarij "rijeka" ──────────────────────────────────────────────────── */

PRS_EROSION_API void prs_erosion_default_river_settings(PrsRiverSettings* out);

/* Zamjenjuje teren dolinom s vijugavim koritom od gole stijene i bazenom na
 * dnu (y = height - 1), briše vodu, sediment i statistiku, i postavlja JEDAN
 * izvor vode na vrhu korita (raniji izvori se uklanjaju). Parametri ostaju.
 * settings == NULL = zadane postavke. Novi teren se čita preko
 * prs_erosion_terrain. Vraća 1 ako uspije. */
PRS_EROSION_API int prs_erosion_setup_river(PrsErosion* sim, const PrsRiverSettings* settings);

/* ── Simulacija ─────────────────────────────────────────────────────────── */

/* `steps` vremenskih koraka zaredom. */
PRS_EROSION_API void prs_erosion_step(PrsErosion* sim, int steps);

/* ── Čitanje stanja ─────────────────────────────────────────────────────── */

PRS_EROSION_API int   prs_erosion_width(const PrsErosion* sim);
PRS_EROSION_API int   prs_erosion_height(const PrsErosion* sim);
PRS_EROSION_API float prs_erosion_cell_size(const PrsErosion* sim);

/* Pokazivači na polja (width*height float-ova, samo za čitanje).
 * VAŽNO: važe samo do sljedećeg poziva prs_erosion_step / reset / destroy,
 * jer simulacija interno zamjenjuje bafere. Uzimati ih svaki frame iznova. */
PRS_EROSION_API const float* prs_erosion_terrain(const PrsErosion* sim);    /* b: visina terena */
PRS_EROSION_API const float* prs_erosion_water(const PrsErosion* sim);      /* d: dubina vode */
PRS_EROSION_API const float* prs_erosion_sediment(const PrsErosion* sim);   /* s: rastvoreni sediment */
PRS_EROSION_API const float* prs_erosion_velocity_x(const PrsErosion* sim); /* u */
PRS_EROSION_API const float* prs_erosion_velocity_y(const PrsErosion* sim); /* v */
/* r: visina stijene; PRS_EROSION_NO_BEDROCK gdje nije zadana. Debljina tla je
 * teren - stijena. Uz rockErodibility > 0 mijenja se tokom simulacije. */
PRS_EROSION_API const float* prs_erosion_bedrock(const PrsErosion* sim);

/* totalWater i totalMaterial se računaju pri pozivu (prolaz preko mreže). */
PRS_EROSION_API void prs_erosion_get_stats(const PrsErosion* sim, PrsErosionStats* out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* PRS_EROSION_H */
