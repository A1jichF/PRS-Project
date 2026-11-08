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

    const std::size_t n = cellCount();
    for (Field* f : {&b, &d, &s, &fL, &fR, &fT, &fB, &u, &v})
        f->assign(n, 0.0f);
}

void Grid::clear()
{
    for (Field* f : {&b, &d, &s, &fL, &fR, &fT, &fB, &u, &v})
        std::fill(f->begin(), f->end(), 0.0f);
}

} // namespace erosion
