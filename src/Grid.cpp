#include "erosion/Grid.h"

#include <algorithm>
#include <stdexcept>

namespace erosion {

Grid::Grid(int w, int h, float cell)
{
    resize(w, h, cell);
}

void Grid::resize(int w, int h, float cell)
{
    if (w <= 0 || h <= 0)
        throw std::invalid_argument("Grid::resize: dimenzije moraju biti pozitivne");
    if (!(cell > 0.0f))
        throw std::invalid_argument("Grid::resize: cellSize mora biti pozitivan");

    width    = w;
    height   = h;
    cellSize = cell;

    // Svako polje na svom pomaku unutar stranice (slotovi 0–9), vidi AlignedAllocator.
    // Stara polja se prvo oslobode, da assign zaista alocira novu memoriju.
    // r je zadnje, da ostala polja zadrže pomake koje su imala prije uvođenja stijene.
    const std::size_t n = cellCount();
    StaggerScope      stagger(0);
    for (Field* f : {&b, &d, &s, &fL, &fR, &fT, &fB, &u, &v}) {
        Field().swap(*f);
        f->assign(n, 0.0f);
    }
    Field().swap(r);
    r.assign(n, kNoBedrock);
}

void Grid::clear()
{
    for (Field* f : {&b, &d, &s, &fL, &fR, &fT, &fB, &u, &v})
        std::fill(f->begin(), f->end(), 0.0f);
    std::fill(r.begin(), r.end(), kNoBedrock);
}

} // namespace erosion
