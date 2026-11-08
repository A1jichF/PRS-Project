#pragma once

#include "erosion/AlignedAllocator.h"

namespace erosion {

struct Grid;

struct FieldRange {
    float min = 0.0f;
    float max = 0.0f;
};

/// Zbir svih vrijednosti polja (akumulacija u double, da greška zaokruživanja
/// ne sakrije stvarni gubitak mase na velikim mrežama).
[[nodiscard]] double fieldSum(const Field& f);

/// Minimum i maksimum polja. Za prazno polje vraća {0, 0}.
[[nodiscard]] FieldRange fieldRange(const Field& f);

/// max |a[i] - b[i]|. Koristimo ga za poređenje sekvencijalne, SIMD i OpenMP verzije.
/// Polja moraju biti iste dužine.
[[nodiscard]] float maxAbsDiff(const Field& a, const Field& b);

/// Ukupna zapremina vode: sum(d) * cellSize².
[[nodiscard]] double totalWater(const Grid& g);

/// Ukupna zapremina materijala: sum(b + s) * cellSize². Mora biti konstantna
/// tokom erozije/depozicije/advekcije.
[[nodiscard]] double totalMaterial(const Grid& g);

/// true ako nijedna vrijednost polja nije NaN ili beskonačna.
[[nodiscard]] bool allFinite(const Field& f);

} // namespace erosion
