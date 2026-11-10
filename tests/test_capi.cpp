// Testovi za C interfejs (prs_erosion.h) kroz pravi prs_erosion.dll.
//
// Glavni test: simulacija kroz DLL mora dati BIT-IDENTIČNE rezultate kao
// statička biblioteka, jer je DLL samo tanak omotač oko istog koda.

#include "TestFramework.h"

#include "prs_erosion.h"

#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/Params.h"
#include "erosion/Simulation.h"
#include "erosion/TerrainGen.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

using namespace erosion;

namespace {

/// max |a[i] - f[i]| za polje iz DLL-a (sirovi pokazivač) i polje iz Grid-a.
float maxDiff(const float* a, const Field& f)
{
    float m = 0.0f;
    for (std::size_t i = 0; i < f.size(); ++i)
        m = std::max(m, std::abs(a[i] - f[i]));
    return m;
}

Grid fbmGrid(int n)
{
    Grid g(n, n, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);
    return g;
}

} // namespace

TEST(capi_verzija_odgovara_headeru)
{
    CHECK(prs_erosion_api_version() == PRS_EROSION_API_VERSION);
}

TEST(capi_zadani_parametri_isti_kao_Params)
{
    PrsErosionParams c;
    std::memset(&c, 0xFF, sizeof c);
    prs_erosion_default_params(&c);
    const Params p;
    CHECK(c.dt == p.dt);
    CHECK(c.gravity == p.gravity);
    CHECK(c.pipeArea == p.pipeArea);
    CHECK(c.rainRate == p.rainRate);
    CHECK(c.evaporationRate == p.evaporationRate);
    CHECK(c.minDepthForVelocity == p.minDepthForVelocity);
    CHECK(c.erosionEnabled == 1);
    CHECK(c.sedimentCapacity == p.sedimentCapacity);
    CHECK(c.dissolvingRate == p.dissolvingRate);
    CHECK(c.depositionRate == p.depositionRate);
    CHECK(c.minTilt == p.minTilt);
    CHECK(c.maxErosionDepth == p.maxErosionDepth);
    CHECK(c.sedimentTransport == PRS_TRANSPORT_FLUX);
}

TEST(capi_dll_bit_identican_statickoj_biblioteci)
{
    // Kiša + izvor vode + erozija, 300 koraka na fBm terenu.
    Grid g = fbmGrid(96);

    Params p;
    p.rainRate = 0.01f;
    Simulation ref(g, p);
    ref.addSource(WaterSource{30, 40, 4.0f, 2.0f});

    PrsErosionParams cp;
    prs_erosion_default_params(&cp);
    cp.rainRate = 0.01f;
    PrsErosion* sim = prs_erosion_create(g.width, g.height, g.cellSize, g.b.data(), &cp);
    CHECK(sim != nullptr);
    if (!sim)
        return;
    CHECK(prs_erosion_add_source(sim, 30, 40, 4.0f, 2.0f) == 1);

    ref.run(300);
    prs_erosion_step(sim, 300);

    const Grid& r = ref.grid();
    CHECK(maxDiff(prs_erosion_terrain(sim), r.b) == 0.0f);
    CHECK(maxDiff(prs_erosion_water(sim), r.d) == 0.0f);
    CHECK(maxDiff(prs_erosion_sediment(sim), r.s) == 0.0f);
    CHECK(maxDiff(prs_erosion_velocity_x(sim), r.u) == 0.0f);
    CHECK(maxDiff(prs_erosion_velocity_y(sim), r.v) == 0.0f);

    PrsErosionStats st;
    prs_erosion_get_stats(sim, &st);
    CHECK(st.steps == 300);
    CHECK(st.simTime == ref.time());
    CHECK(st.totalWater == totalWater(r));
    CHECK(st.totalMaterial == totalMaterial(r));
    CHECK(st.lastStepMs > 0.0);

    prs_erosion_destroy(sim);
}

TEST(capi_promjena_parametara_usred_simulacije)
{
    // Isto kao kad renderer mijenja slider: set_params između koraka.
    Grid g = fbmGrid(64);
    Params p;
    p.rainRate = 0.02f;
    Simulation ref(g, p);

    PrsErosionParams cp;
    prs_erosion_default_params(&cp);
    cp.rainRate = 0.02f;
    PrsErosion* sim = prs_erosion_create(g.width, g.height, g.cellSize, g.b.data(), &cp);
    CHECK(sim != nullptr);
    if (!sim)
        return;

    ref.run(100);
    prs_erosion_step(sim, 100);

    ref.params().sedimentTransport = SedimentTransport::SemiLagrangian;
    ref.params().depositionRate    = 0.3f;
    cp.sedimentTransport = PRS_TRANSPORT_SEMI_LAGRANGIAN;
    cp.depositionRate    = 0.3f;
    CHECK(prs_erosion_set_params(sim, &cp) == 1);

    PrsErosionParams back;
    prs_erosion_get_params(sim, &back);
    CHECK(back.sedimentTransport == PRS_TRANSPORT_SEMI_LAGRANGIAN);
    CHECK(back.depositionRate == 0.3f);

    ref.run(100);
    prs_erosion_step(sim, 100);
    CHECK(maxDiff(prs_erosion_terrain(sim), ref.grid().b) == 0.0f);
    CHECK(maxDiff(prs_erosion_sediment(sim), ref.grid().s) == 0.0f);

    prs_erosion_destroy(sim);
}

TEST(capi_reset_vraca_pocetno_stanje)
{
    Grid g = fbmGrid(48);
    PrsErosionParams cp;
    prs_erosion_default_params(&cp);
    cp.rainRate = 0.05f;
    PrsErosion* sim = prs_erosion_create(g.width, g.height, g.cellSize, g.b.data(), &cp);
    CHECK(sim != nullptr);
    if (!sim)
        return;
    CHECK(prs_erosion_add_source(sim, 10, 10, 2.0f, 1.0f) == 1);

    prs_erosion_step(sim, 50);
    const std::vector<float> after50(prs_erosion_terrain(sim),
                                     prs_erosion_terrain(sim) + g.cellCount());

    CHECK(prs_erosion_reset(sim, g.b.data()) == 1);
    PrsErosionStats st;
    prs_erosion_get_stats(sim, &st);
    CHECK(st.steps == 0);
    CHECK(st.totalWater == 0.0);
    CHECK(maxDiff(prs_erosion_terrain(sim), g.b) == 0.0f);

    // Parametri i izvor su zadržani, pa ponovljenih 50 koraka daje isto.
    prs_erosion_step(sim, 50);
    float m = 0.0f;
    for (std::size_t i = 0; i < after50.size(); ++i)
        m = std::max(m, std::abs(after50[i] - prs_erosion_terrain(sim)[i]));
    CHECK(m == 0.0f);

    prs_erosion_destroy(sim);
}

TEST(capi_greske_ne_ruse_program)
{
    CHECK(prs_erosion_create(1, 10, 1.0f, nullptr, nullptr) == nullptr);
    CHECK(std::strlen(prs_erosion_last_error()) > 0);
    CHECK(prs_erosion_create(10, 10, 0.0f, nullptr, nullptr) == nullptr);

    PrsErosionParams bad;
    prs_erosion_default_params(&bad);
    bad.dt = 0.0f;
    CHECK(prs_erosion_create(10, 10, 1.0f, nullptr, &bad) == nullptr);
    bad.dt                = 0.02f;
    bad.sedimentTransport = 7;
    CHECK(prs_erosion_create(10, 10, 1.0f, nullptr, &bad) == nullptr);

    // Neispravni parametri se odbiju i stari ostaju.
    PrsErosion* sim = prs_erosion_create(10, 10, 1.0f, nullptr, nullptr);
    CHECK(sim != nullptr);
    CHECK(std::strlen(prs_erosion_last_error()) == 0);
    bad.dt = -1.0f;
    CHECK(prs_erosion_set_params(sim, &bad) == 0);
    PrsErosionParams cur;
    prs_erosion_get_params(sim, &cur);
    CHECK(cur.dt == Params{}.dt);
    CHECK(cur.sedimentTransport == PRS_TRANSPORT_FLUX);

    // NULL handle: ništa se ne desi.
    prs_erosion_step(nullptr, 10);
    prs_erosion_destroy(nullptr);
    CHECK(prs_erosion_terrain(nullptr) == nullptr);
    CHECK(prs_erosion_width(nullptr) == 0);
    CHECK(prs_erosion_reset(nullptr, nullptr) == 0);

    CHECK(prs_erosion_width(sim) == 10);
    CHECK(prs_erosion_height(sim) == 10);
    CHECK(prs_erosion_cell_size(sim) == 1.0f);
    prs_erosion_destroy(sim);
}
