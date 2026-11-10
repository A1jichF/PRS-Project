#pragma once

// Početna stanja za testove i demonstracije.

#include "erosion/Kernels.h"

namespace erosion {

struct Grid;

/// Ravan teren visine `terrainHeight` prekriven vodom dubine `waterDepth`.
/// Brišu se voda, sediment, fluks i brzina. Test: voda mora ostati mirna.
void setupFlatLake(Grid& g, float terrainHeight, float waterDepth);

/// Napuni teren vodom do nivoa `level`: d = max(0, level - b).
/// Površina vode je ravna, pa i na neravnom terenu voda mora ostati mirna.
/// Teren (b) se ne mijenja; fluks i brzina se brišu.
void fillToLevel(Grid& g, float level);

/// Dam-break (pucanje brane): ravan teren na visini 0, voda dubine `baseDepth`
/// svuda i stub vode dubine `columnDepth` u krugu poluprečnika `radius`
/// (u ćelijama) oko centra mreže. Stub je simetričan po x, y i dijagonali.
void setupDamBreak(Grid& g, float baseDepth, float columnDepth, int radius);

/// Parametri scenarija "rijeka" (setupRiverValley). Dužine i visine su u
/// svjetskim jedinicama, osim gdje piše da su u ćelijama ili udjelima mreže.
struct RiverSettings {
    float slope            = 0.12f; ///< Pad doline nizvodno (visina po jedinici dužine)
    float valleySlope      = 0.08f; ///< Nagib strana doline prema koritu
    float channelHalfWidth = 4.0f;  ///< Poluširina ravnog dna korita [ćelije]
    float channelDepth     = 3.0f;  ///< Dubina korita ispod obale
    float bankWidth        = 2.0f;  ///< Širina obale (zida korita) [ćelije]
    float meanderAmplitude = 0.10f; ///< Amplituda vijuganja kao udio širine mreže
    float meanders         = 1.25f; ///< Broj punih vijuga od izvora do bazena
    float basinFraction    = 0.22f; ///< Donji dio mreže (udio visine) je bazen bez korita
    float soilDepth        = 1.5f;  ///< Debljina tla izvan korita (u koritu je gola stijena)
    float sourceRate       = 0.6f;  ///< Izvor: visina vode u sekundi po ćeliji izvora
};

/// Scenarij "rijeka": dolina koja pada od y = 0 prema y = height - 1, s
/// vijugavim koritom usječenim u stijenu i bazenom na dnu u kojem se voda
/// skuplja u jezero (rubovi mreže su zatvoreni).
///
/// Postavlja teren b i stijenu r: dno i obale korita su gola stijena (r = b),
/// pa se korito ne može produbiti ni proširiti; izvan korita je sloj tla
/// debljine RiverSettings::soilDepth. Voda, sediment, fluks i brzina se brišu.
///
/// Vraća izvor vode na vrhu korita; treba ga dodati simulaciji
/// (Simulation::addSource), inače rijeka nema vode.
[[nodiscard]] WaterSource setupRiverValley(Grid& g, const RiverSettings& settings = {});

/// Maska korita iz setupRiverValley za iste postavke: 1 na dnu i obalama
/// korita, 0 drugdje (g.cellCount() vrijednosti). Za testove i prikaz.
void riverChannelMask(const Grid& g, const RiverSettings& settings, float* mask);

} // namespace erosion
