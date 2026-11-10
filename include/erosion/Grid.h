#pragma once

#include "erosion/AlignedAllocator.h"

#include <cstddef>
#include <limits>

namespace erosion {

/// Vrijednost polja Grid::r koja znači "ovdje nema stijene": najmanji konačni
/// float, pa je debljina tla b - r praktično neograničena i erozija se ponaša
/// kao da stijena ne postoji. Namjerno nije -inf, da polje ostane konačno
/// (allFinite, a kasnije i SIMD s -ffinite-math-only).
inline constexpr float kNoBedrock = std::numeric_limits<float>::lowest();

/// Stanje simulacije na mreži width × height ćelija, u SoA rasporedu.
///
/// Konvencije (iste kao u opengl-forward-renderer, radi kasnije integracije):
///   - indeks ćelije: idx(x, y) = y * width + x (red po red, row-major)
///   - x raste udesno, y raste "naniže" (red 0 je prvi red u PGM slici);
///     u 3D svijetu renderera y mreže odgovara osi Z, a visine idu po osi Y
///   - sve veličine su float, visine u svjetskim jedinicama
///
/// Susjedi za fluks:  L = (x-1, y),  R = (x+1, y),  T = (x, y-1),  B = (x, y+1).
///
/// Pomoćne bafere za double-buffering (npr. novi sediment pri advekciji)
/// drži simulacija, ne Grid: Grid opisuje samo trenutno stanje.
struct Grid {
    int   width    = 0;
    int   height   = 0;
    float cellSize = 1.0f; ///< Razmak između centara ćelija (dužina cijevi l)

    Field b;  ///< Visina terena
    Field d;  ///< Visina (dubina) vode
    Field s;  ///< Količina rastvorenog sedimenta

    Field fL; ///< Odlazni fluks prema lijevom susjedu
    Field fR; ///< Odlazni fluks prema desnom susjedu
    Field fT; ///< Odlazni fluks prema gornjem susjedu (y-1)
    Field fB; ///< Odlazni fluks prema donjem susjedu (y+1)

    Field u;  ///< Brzina vode po x
    Field v;  ///< Brzina vode po y

    /// Visina vrha stijene (bedrock). Između r i b je sloj tla debljine b - r,
    /// koji voda može odnijeti; stijenu ne može (ili može sporije, zavisno od
    /// tvrdoće). Gdje je r == b, stijena je gola. Zadano kNoBedrock: nema stijene.
    Field r;

    Grid() = default;
    Grid(int width, int height, float cellSize);

    /// Alocira sva polja na width × height i postavlja ih na nulu
    /// (stijenu r na kNoBedrock).
    void resize(int width, int height, float cellSize);

    /// Postavlja sva polja na nulu, a stijenu r na kNoBedrock (dimenzije ostaju iste).
    void clear();

    [[nodiscard]] std::size_t cellCount() const
    {
        return static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    }

    [[nodiscard]] std::size_t idx(int x, int y) const
    {
        return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
               static_cast<std::size_t>(x);
    }

    [[nodiscard]] bool inBounds(int x, int y) const
    {
        return x >= 0 && y >= 0 && x < width && y < height;
    }
};

} // namespace erosion
