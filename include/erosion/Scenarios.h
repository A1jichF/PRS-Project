#pragma once

// Početna stanja za testove i demonstracije.

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

} // namespace erosion
