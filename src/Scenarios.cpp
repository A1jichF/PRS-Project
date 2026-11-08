#include "erosion/Scenarios.h"

#include "erosion/Grid.h"
#include "erosion/TerrainGen.h"

#include <algorithm>

namespace erosion {

namespace {

void clearDynamics(Grid& g)
{
    for (Field* f : {&g.s, &g.fL, &g.fR, &g.fT, &g.fB, &g.u, &g.v})
        fillConstant(*f, 0.0f);
}

} // namespace

void setupFlatLake(Grid& g, float terrainHeight, float waterDepth)
{
    fillConstant(g.b, terrainHeight);
    fillConstant(g.d, waterDepth);
    clearDynamics(g);
}

void fillToLevel(Grid& g, float level)
{
    for (std::size_t i = 0; i < g.cellCount(); ++i)
        g.d[i] = std::max(0.0f, level - g.b[i]);
    clearDynamics(g);
}

void setupDamBreak(Grid& g, float baseDepth, float columnDepth, int radius)
{
    setupFlatLake(g, 0.0f, baseDepth);

    // Centar kao dvostruka koordinata, da bude tačan i za parne dimenzije.
    const int cx2 = g.width - 1;
    const int cy2 = g.height - 1;
    const int r2  = 4 * radius * radius;
    for (int y = 0; y < g.height; ++y) {
        for (int x = 0; x < g.width; ++x) {
            const int dx = 2 * x - cx2;
            const int dy = 2 * y - cy2;
            if (dx * dx + dy * dy <= r2)
                g.d[g.idx(x, y)] = columnDepth;
        }
    }
}

} // namespace erosion
