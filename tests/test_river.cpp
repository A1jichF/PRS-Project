// Testovi za scenarij "rijeka": dolina s koritom od stijene (setupRiverValley).

#include "TestFramework.h"

#include "erosion/Bedrock.h"
#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/Params.h"
#include "erosion/Scenarios.h"
#include "erosion/Simulation.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace erosion;

namespace {

constexpr int kSize  = 128;
constexpr int kSteps = 4000;

struct RiverRun {
    Field              b0;      ///< Početni teren
    std::vector<float> channel; ///< 1 u koritu, 0 drugdje
    double             material0 = 0.0;
    Simulation         sim;
};

RiverRun runRiver(bool withBedrock, int steps)
{
    Grid                g(kSize, kSize, 2.0f);
    const RiverSettings settings;
    const WaterSource   source = setupRiverValley(g, settings);
    std::vector<float>  channel(g.cellCount());
    riverChannelMask(g, settings, channel.data());
    if (!withBedrock)
        clearBedrock(g);

    RiverRun run{g.b, std::move(channel), totalMaterial(g), Simulation(std::move(g))};
    run.sim.addSource(source);
    run.sim.run(steps);
    return run;
}

} // namespace

TEST(rijeka_korito_je_gola_stijena_a_okolo_je_tlo)
{
    Grid                g(kSize, kSize, 2.0f);
    const RiverSettings settings;
    const WaterSource   source = setupRiverValley(g, settings);
    std::vector<float>  channel(g.cellCount());
    riverChannelMask(g, settings, channel.data());

    std::size_t channelCells = 0;
    for (std::size_t i = 0; i < g.cellCount(); ++i) {
        if (channel[i] > 0.5f) {
            ++channelCells;
            CHECK(g.r[i] == g.b[i]);
        } else {
            CHECK(g.r[i] == g.b[i] - settings.soilDepth);
        }
    }
    // Korito postoji, ali je mali dio mreže.
    CHECK(channelCells > 100);
    CHECK(channelCells < g.cellCount() / 5);

    // Izvor je u koritu, a početno stanje je suho i mirno.
    CHECK(g.inBounds(source.x, source.y));
    CHECK(channel[g.idx(source.x, source.y)] > 0.5f);
    CHECK(source.rate == settings.sourceRate);
    CHECK(fieldSum(g.d) == 0.0);
    CHECK(fieldSum(g.s) == 0.0);
    CHECK(allFinite(g.b));
    CHECK(fieldRange(g.b).min >= 0.0f);
}

TEST(rijeka_dno_korita_pada_nizvodno)
{
    // Najniža tačka svakog reda (dno korita) ne smije rasti nizvodno, inače
    // bi voda stajala u koritu umjesto da teče u bazen.
    Grid g(kSize, kSize, 2.0f);
    (void)setupRiverValley(g, RiverSettings{});

    float previous = fieldRange(g.b).max;
    for (int y = 0; y < g.height; ++y) {
        float rowMin = g.b[g.idx(0, y)];
        for (int x = 1; x < g.width; ++x)
            rowMin = std::min(rowMin, g.b[g.idx(x, y)]);
        CHECK(rowMin <= previous);
        previous = rowMin;
    }
    CHECK(previous == 0.0f); // bazen je na visini 0
}

TEST(rijeka_korito_od_stijene_se_ne_produbljuje)
{
    const RiverRun run = runRiver(/*withBedrock=*/true, kSteps);
    const Grid&    g   = run.sim.grid();

    // U koritu teren nikad nije niži od početnog (može biti viši: talog).
    float deepest = 0.0f, fastest = 0.0f;
    for (std::size_t i = 0; i < g.cellCount(); ++i) {
        if (run.channel[i] < 0.5f)
            continue;
        deepest = std::min(deepest, g.b[i] - run.b0[i]);
        fastest = std::max(fastest, std::hypot(g.u[i], g.v[i]));
    }
    CHECK(deepest == 0.0f);
    CHECK(fastest > 1.0f); // rijeka teče

    // Voda je stigla do bazena na dnu i tu se skuplja (rubovi su zatvoreni).
    double basinWater = 0.0;
    for (int y = g.height - g.height / 6; y < g.height; ++y)
        for (int x = 0; x < g.width; ++x)
            basinWater += g.d[g.idx(x, y)];
    CHECK(basinWater > 0.0);

    CHECK(std::abs(totalMaterial(g) - run.material0) / run.material0 < 1e-5);
    CHECK(allFinite(g.b));
    CHECK(allFinite(g.d));
}

TEST(rijeka_bez_stijene_se_usijeca)
{
    // Ista rijeka bez stijene: korito se produbljuje. Ovo pokazuje da test
    // iznad zaista provjerava stijenu, a ne preslabu eroziju.
    const RiverRun run = runRiver(/*withBedrock=*/false, kSteps);
    const Grid&    g   = run.sim.grid();

    float deepest = 0.0f;
    for (std::size_t i = 0; i < g.cellCount(); ++i)
        if (run.channel[i] > 0.5f)
            deepest = std::min(deepest, g.b[i] - run.b0[i]);
    CHECK(deepest < -1.0f);
}
