#pragma once

namespace erosion {

/// Kako sediment putuje s vodom (korak 5).
enum class SedimentTransport {
    /// Konzervativno, kroz iste virtualne cijevi kao voda: ćelija šalje susjedu
    /// isti udio svog sedimenta kao i vode. Masa se čuva tačno.
    /// Odstupanje od Mei et al. (vidi README, korak 5).
    Flux,
    /// Mei et al. 2007: semi-Lagrangian advekcija s(x - v*dt). Ne čuva masu:
    /// gdje se tok usporava i skuplja (dno padine), sediment nestaje.
    SemiLagrangian,
};

/// Parametri simulacije (Mei, Decaudin, Hu 2007: virtual pipe model).
///
/// Razmak ćelija (dužina cijevi l) nije ovdje nego u Grid::cellSize, da ne
/// postoji na dva mjesta. Konstante erozije (Kc, Ks, Kd, Ke) su početne
/// vrijednosti i podešavat ćemo ih u koracima 4–6.
struct Params {
    float dt = 0.02f; ///< Vremenski korak [s]

    // ── Virtualne cijevi ─────────────────────────────────────────────────────
    float gravity  = 9.81f; ///< g
    /// A, poprečni presjek cijevi. Napomena: brzina talasa u pipe modelu ne
    /// zavisi od stvarne dubine vode, nego se ponaša kao da je dubina A / l.
    float pipeArea = 1.0f;

    // ── Voda ─────────────────────────────────────────────────────────────────
    float rainRate        = 0.0f;   ///< Visina vode dodana svakoj ćeliji u sekundi
    float evaporationRate = 0.015f; ///< Ke: udio vode koji ispari u sekundi (korak 6)
    /// Ispod ove prosječne dubine brzina se postavlja na 0, jer bi dijeljenje
    /// malom dubinom dalo ogromne (nefizikalne) brzine na skoro suhim ćelijama.
    float minDepthForVelocity = 1e-4f;

    // ── Erozija i depozicija (korak 4) ───────────────────────────────────────
    //   C = Kc * max(sin α, minTilt) * |v| * lmax(d)
    //   C > s: erozija    b -= Ks*dt*(C - s),  s += isto
    //   C ≤ s: depozicija b += Kd*dt*(s - C),  s -= isto
    //
    // Razlika u odnosu na Mei et al.: tamo su Ks i Kd udio po koraku, a ovdje
    // su brzine [1/s] pomnožene s dt (ograničeno na 1), da rezultat ne zavisi
    // od izbora dt.
    bool  erosionEnabled   = true;
    float sedimentCapacity = 0.5f;  ///< Kc
    float dissolvingRate   = 0.5f;  ///< Ks [1/s]
    float depositionRate   = 1.0f;  ///< Kd [1/s]
    float minTilt          = 0.05f; ///< Donja granica za sin α, da C ne padne na 0 na ravnom
    /// Kdmax iz Šťava et al. 2008: lmax(d) raste linearno od 0 (d = 0) do 1
    /// (d >= Kdmax). Tanak sloj vode tako ne erodira punom snagom, što bi
    /// inače radio zbog umjetno velike brzine na ćelijama koje se isprazne
    /// u jednom koraku (u = l / dt). 0 isključuje faktor.
    float maxErosionDepth  = 0.1f;

    // ── Transport sedimenta (korak 5) ────────────────────────────────────────
    SedimentTransport sedimentTransport = SedimentTransport::Flux;

    // ── Termalna erozija (opcionalno) ────────────────────────────────────────
    float talusAngleDeg = 38.0f; ///< Ugao talusa
    float thermalRate   = 0.15f; ///< Udio viška materijala koji se prebaci po koraku
};

/// Parametri za pravljenje terena (ne utiču na samu simulaciju).
/// Zadane vrijednosti odgovaraju opengl-forward-rendereru: 512 jedinica
/// širine na 256 ćelija (cellSize = 2) i maksimalna visina 80.
struct WorldSettings {
    float cellSize    = 2.0f;
    float heightScale = 80.0f;
};

} // namespace erosion
