#include "erosion/Bedrock.h"

#include "erosion/Grid.h"

#include <algorithm>
#include <cmath>

namespace erosion {

void clearBedrock(Grid& g)
{
    std::fill(g.r.begin(), g.r.end(), kNoBedrock);
}

bool hasBedrock(const Grid& g)
{
    return std::any_of(g.r.begin(), g.r.end(), [](float r) { return r != kNoBedrock; });
}

void setBedrockUniformSoil(Grid& g, float soilDepth)
{
    if (!std::isfinite(soilDepth)) { // beskonačno tlo = nema stijene
        clearBedrock(g);
        return;
    }
    const float       depth = std::max(0.0f, soilDepth);
    const std::size_t n     = g.cellCount();
    for (std::size_t i = 0; i < n; ++i)
        g.r[i] = g.b[i] - depth;
}

void setBedrockFromSoilDepth(Grid& g, const float* soilDepth)
{
    const std::size_t n = g.cellCount();
    for (std::size_t i = 0; i < n; ++i) {
        const float depth = soilDepth[i];
        g.r[i] = std::isfinite(depth) ? g.b[i] - std::max(0.0f, depth) : kNoBedrock;
    }
}

void setBedrockLevel(Grid& g, float level)
{
    if (!(level > kNoBedrock)) { // NaN ili -inf: nema stijene
        clearBedrock(g);
        return;
    }
    const std::size_t n = g.cellCount();
    for (std::size_t i = 0; i < n; ++i)
        g.r[i] = std::min(g.b[i], level);
}

void setBedrockHeights(Grid& g, const float* heights)
{
    const std::size_t n = g.cellCount();
    for (std::size_t i = 0; i < n; ++i) {
        const float h = heights[i];
        g.r[i] = std::isfinite(h) ? std::min(g.b[i], h) : kNoBedrock;
    }
}

BedrockStats bedrockStats(const Grid& g, float bareThreshold)
{
    BedrockStats      stats;
    const std::size_t n = g.cellCount();
    double            soil = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        if (g.r[i] == kNoBedrock)
            continue;
        const float depth = std::max(0.0f, g.b[i] - g.r[i]);
        ++stats.rockCells;
        if (depth <= bareThreshold)
            ++stats.bareCells;
        soil += static_cast<double>(depth);
    }
    stats.soilVolume = soil * static_cast<double>(g.cellSize) * static_cast<double>(g.cellSize);
    return stats;
}

} // namespace erosion
