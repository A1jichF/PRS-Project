// Testovi za korak 6: isparavanje.

#include "TestFramework.h"

#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/Kernels.h"
#include "erosion/Params.h"
#include "erosion/Scenarios.h"
#include "erosion/Simulation.h"
#include "erosion/TerrainGen.h"

#include <cmath>

using namespace erosion;

namespace {

double relErr(double actual, double expected)
{
    return std::abs(actual - expected) / std::abs(expected);
}

} // namespace

TEST(isparavanje_eksponencijalni_pad)
{
    // Ravno jezero isparava jednako svuda, pa ostaje mirno: d = d0 * (1 - Ke*dt)^N.
    Grid g(32, 32, 2.0f);
    setupFlatLake(g, 5.0f, 2.0f);
    Params p;
    p.evaporationRate = 0.5f;
    Simulation sim(std::move(g), p);
    sim.run(100);

    const double expected = 2.0 * std::pow(1.0 - 0.5 * p.dt, 100);
    CHECK(relErr(fieldRange(sim.grid().d).min, expected) < 1e-5);
    CHECK(relErr(fieldRange(sim.grid().d).max, expected) < 1e-5);
}

TEST(ukupna_voda_pada_tacno_faktorom_i_kad_voda_tece)
{
    // Tok ne mijenja ukupnu vodu, a isparavanje je množi istim faktorom svuda,
    // pa ukupna voda nakon N koraka mora biti W0 * f^N bez obzira na tok.
    Grid g(64, 64, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);
    fillConstant(g.d, 1.0f);
    const double w0 = totalWater(g);

    Params p;
    p.evaporationRate = 0.1f;
    p.erosionEnabled  = false;
    Simulation sim(std::move(g), p);
    sim.run(300);

    const double expected = w0 * std::pow(1.0 - 0.1 * p.dt, 300);
    CHECK(relErr(totalWater(sim.grid()), expected) < 1e-4);
}

TEST(preveliko_isparavanje_ne_daje_negativnu_vodu)
{
    Grid g(8, 8, 1.0f);
    fillConstant(g.d, 1.0f);
    Params p;
    p.evaporationRate = 1000.0f; // Ke * dt = 20 > 1
    seq::evaporate(g, p);
    CHECK(fieldRange(g.d).min == 0.0f);
    CHECK(fieldRange(g.d).max == 0.0f);
}

TEST(isparavanje_ne_dira_teren_ni_sediment)
{
    Grid g(8, 8, 1.0f);
    fillConstant(g.b, 3.0f);
    fillConstant(g.d, 1.0f);
    fillConstant(g.s, 0.4f);
    Params p;
    p.evaporationRate = 5.0f;
    seq::evaporate(g, p);
    CHECK(fieldRange(g.b).min == 3.0f && fieldRange(g.b).max == 3.0f);
    CHECK(fieldRange(g.s).min == 0.4f && fieldRange(g.s).max == 0.4f);
    CHECK(fieldRange(g.d).max < 1.0f);
}

TEST(kisa_i_isparavanje_dostignu_ravnotezu)
{
    // Svaki korak: d = (d + r*dt) * f, f = 1 - Ke*dt. Ravnoteža: d* = r*dt*f / (1 - f).
    Grid g(16, 16, 1.0f);
    Params p;
    p.rainRate        = 0.1f;
    p.evaporationRate = 0.5f;
    Simulation sim(std::move(g), p);
    sim.run(3000);

    const double f        = 1.0 - 0.5 * p.dt;
    const double expected = 0.1 * p.dt * f / (1.0 - f); // 0.198
    CHECK(relErr(fieldRange(sim.grid().d).min, expected) < 1e-4);
    CHECK(relErr(fieldRange(sim.grid().d).max, expected) < 1e-4);
}

TEST(materijal_se_cuva_uz_isparavanje)
{
    Grid g(128, 128, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);
    const double m0 = totalMaterial(g);

    Params p; // zadano isparavanje 0.015
    p.rainRate = 0.02f;
    Simulation sim(std::move(g), p);
    sim.run(3000);

    const Grid& r = sim.grid();
    CHECK(relErr(totalMaterial(r), m0) < 1e-5);
    CHECK(allFinite(r.b) && allFinite(r.s) && allFinite(r.d));
    CHECK(fieldRange(r.s).min >= 0.0f);
    CHECK(fieldRange(r.d).min >= 0.0f);
}
