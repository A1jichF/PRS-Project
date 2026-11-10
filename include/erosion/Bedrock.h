#pragma once

// Zadavanje stijene (bedrock, polje Grid::r) i statistika sloja tla.
//
// Stijena se zadaje POSLIJE terena (b), jer se sve funkcije oslanjaju na
// trenutnu visinu terena i garantuju r <= b. Ako se teren kasnije promijeni
// (novi fBm, učitan PGM), stijenu treba zadati ponovo.

#include <cstddef>

namespace erosion {

struct Grid;

/// Uklanja stijenu: r = kNoBedrock svuda (erozija opet nije ograničena).
void clearBedrock(Grid& g);

/// true ako je stijena zadana u barem jednoj ćeliji.
[[nodiscard]] bool hasBedrock(const Grid& g);

/// Isti sloj tla svuda: r = b - soilDepth. Stijena tada prati oblik terena.
/// soilDepth = 0 daje golu stijenu svuda; negativna vrijednost se računa kao 0,
/// a NaN ili beskonačna uklanja stijenu.
void setBedrockUniformSoil(Grid& g, float soilDepth);

/// Sloj tla zadan po ćeliji: r[i] = b[i] - soilDepth[i]. `soilDepth` ima
/// g.cellCount() vrijednosti. Negativna debljina se računa kao 0 (gola
/// stijena), a NaN ili beskonačna kao "ovdje nema stijene".
void setBedrockFromSoilDepth(Grid& g, const float* soilDepth);

/// Stijena do visine `level`: r = min(b, level). Ćelije na visini <= level
/// postaju gola stijena, a iznad toga je tlo debljine b - level.
/// NaN ili -inf uklanja stijenu.
void setBedrockLevel(Grid& g, float level);

/// Proizvoljna mapa visina stijene (g.cellCount() vrijednosti, npr. iz PGM-a
/// ili iz renderera). Stijena viša od terena se spušta na teren (r = b), a
/// NaN ili beskonačna vrijednost znači "ovdje nema stijene".
void setBedrockHeights(Grid& g, const float* heights);

struct BedrockStats {
    std::size_t rockCells  = 0;   ///< Ćelije u kojima je stijena zadana
    std::size_t bareCells  = 0;   ///< Od toga gola stijena (tlo tanje od bareThreshold)
    double      soilVolume = 0.0; ///< Zapremina tla iznad stijene: sum(b - r) * cellSize²
};

/// Statistika sloja tla. Ćelije bez stijene (kNoBedrock) se ne broje, jer je
/// njihovo tlo po definiciji neograničeno.
[[nodiscard]] BedrockStats bedrockStats(const Grid& g, float bareThreshold = 0.0f);

} // namespace erosion
