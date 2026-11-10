#include "erosion/Scenarios.h"

#include "erosion/Grid.h"
#include "erosion/TerrainGen.h"

#include <algorithm>
#include <cmath>

namespace erosion {

namespace {

void clearDynamics(Grid& g)
{
    for (Field* f : {&g.s, &g.fL, &g.fR, &g.fT, &g.fB, &g.u, &g.v})
        fillConstant(*f, 0.0f);
}

float smoothstep(float edge0, float edge1, float x)
{
    const float t = std::clamp((x - edge0) / (edge1 - edge0), 0.0f, 1.0f);
    return t * t * (3.0f - 2.0f * t);
}

/// Geometrija doline iz setupRiverValley, zajednička za teren i masku korita.
struct RiverGeometry {
    RiverGeometry(const Grid& g, const RiverSettings& s)
        : settings(s)
        , midX(0.5f * static_cast<float>(g.width - 1))
        , amplitude(s.meanderAmplitude * static_cast<float>(g.width - 1))
        , basinRow(std::max(1.0f, (1.0f - std::clamp(s.basinFraction, 0.0f, 0.9f)) *
                                      static_cast<float>(g.height - 1)))
        , fadeRows(std::max(1.0f, 0.08f * static_cast<float>(g.height - 1)))
    {
    }

    /// x koordinata sredine korita u redu y [ćelije]. U bazenu ostaje gdje je stala.
    [[nodiscard]] float centerX(float y) const
    {
        constexpr float kTwoPi = 6.283185307179586f;
        const float     t      = std::clamp(y / basinRow, 0.0f, 1.0f);
        return midX + amplitude * std::sin(kTwoPi * settings.meanders * t);
    }

    /// 1 uzvodno, glatko pada na 0 pred bazenom: korito se tu postepeno gubi.
    [[nodiscard]] float channelFade(float y) const
    {
        return 1.0f - smoothstep(basinRow - fadeRows, basinRow, y);
    }

    /// Udaljenost ćelije od ruba ravnog dna korita [ćelije]; 0 na dnu.
    ///
    /// Vodoravna udaljenost od sredine se množi kosinusom ugla korita prema osi
    /// y, što daje (približno) udaljenost okomito na korito. Bez toga bi korito
    /// u krivinama, gdje ide ukoso, bilo uže nego na ravnim dijelovima.
    [[nodiscard]] float outsideBed(float x, float y) const
    {
        const float dxdy = centerX(y + 0.5f) - centerX(y - 0.5f);
        const float cosA = 1.0f / std::sqrt(1.0f + dxdy * dxdy);
        return std::max(0.0f, std::abs(x - centerX(y)) * cosA - settings.channelHalfWidth);
    }

    /// true na dnu i obalama korita (tu je stijena gola).
    [[nodiscard]] bool inChannel(float x, float y) const
    {
        return channelFade(y) > 0.0f && outsideBed(x, y) <= settings.bankWidth;
    }

    RiverSettings settings;
    float         midX;
    float         amplitude;
    float         basinRow; ///< Red u kojem počinje bazen
    float         fadeRows; ///< Dužina prelaza korito -> bazen
};

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

WaterSource setupRiverValley(Grid& g, const RiverSettings& settings)
{
    const RiverGeometry geo(g, settings);
    const float         l         = g.cellSize;
    const float         bankWidth = std::max(1e-3f, settings.bankWidth);
    const float         soilDepth = std::max(0.0f, settings.soilDepth);

    for (int y = 0; y < g.height; ++y) {
        const float fy    = static_cast<float>(y);
        const float fade  = geo.channelFade(fy);
        // Dno korita pada nizvodno do bazena, a bazen je na visini baseHeight.
        const float along = settings.baseHeight +
                            settings.slope * l * std::max(0.0f, geo.basinRow - fy);
        // Strane bazena su blaže od strana doline, da se voda razlije u jezero.
        const float sideScale = 0.3f + 0.7f * fade;

        for (int x = 0; x < g.width; ++x) {
            const float outside = geo.outsideBed(static_cast<float>(x), fy);
            // Obala: od dna korita (0) do pune dubine korita preko bankWidth ćelija.
            const float bank = settings.channelDepth * fade * smoothstep(0.0f, bankWidth, outside);
            const float side = settings.valleySlope * sideScale * l *
                               std::max(0.0f, outside - settings.bankWidth);

            const std::size_t i = g.idx(x, y);
            g.b[i] = along + bank + side;
            g.r[i] = geo.inChannel(static_cast<float>(x), fy) ? g.b[i] : g.b[i] - soilDepth;
        }
    }

    fillConstant(g.d, 0.0f);
    clearDynamics(g);

    // Izvor na vrhu korita, malo odmaknut od gornjeg ruba mreže.
    const int sourceY = std::min(g.height - 1,
                                 static_cast<int>(std::ceil(settings.channelHalfWidth)) + 2);
    WaterSource source;
    source.x      = std::clamp(static_cast<int>(std::lround(geo.centerX(static_cast<float>(sourceY)))),
                               0, g.width - 1);
    source.y      = sourceY;
    source.radius = std::max(1.0f, settings.channelHalfWidth);
    source.rate   = settings.sourceRate;
    return source;
}

void riverChannelMask(const Grid& g, const RiverSettings& settings, float* mask)
{
    const RiverGeometry geo(g, settings);
    for (int y = 0; y < g.height; ++y)
        for (int x = 0; x < g.width; ++x)
            mask[g.idx(x, y)] =
                geo.inChannel(static_cast<float>(x), static_cast<float>(y)) ? 1.0f : 0.0f;
}

} // namespace erosion
