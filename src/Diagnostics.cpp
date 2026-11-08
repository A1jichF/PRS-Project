#include "erosion/Diagnostics.h"

#include "erosion/Grid.h"
#include "erosion/Params.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace erosion {

double fieldSum(const Field& f)
{
    double sum = 0.0;
    for (float x : f)
        sum += static_cast<double>(x);
    return sum;
}

FieldRange fieldRange(const Field& f)
{
    if (f.empty())
        return {};
    const auto [mn, mx] = std::minmax_element(f.begin(), f.end());
    return {*mn, *mx};
}

float maxAbsDiff(const Field& a, const Field& b)
{
    if (a.size() != b.size())
        throw std::invalid_argument("maxAbsDiff: polja nisu iste dužine");
    float m = 0.0f;
    for (std::size_t i = 0; i < a.size(); ++i)
        m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

double totalWater(const Grid& g)
{
    const double area = static_cast<double>(g.cellSize) * g.cellSize;
    return fieldSum(g.d) * area;
}

double totalMaterial(const Grid& g)
{
    const double area = static_cast<double>(g.cellSize) * g.cellSize;
    return (fieldSum(g.b) + fieldSum(g.s)) * area;
}

bool allFinite(const Field& f)
{
    return std::all_of(f.begin(), f.end(), [](float x) { return std::isfinite(x); });
}

CflInfo velocityCfl(const Grid& g, const Params& p)
{
    float m = 0.0f;
    for (std::size_t i = 0; i < g.cellCount(); ++i)
        m = std::max(m, std::max(std::abs(g.u[i]), std::abs(g.v[i])));
    return {m, p.dt * m / g.cellSize};
}

} // namespace erosion
