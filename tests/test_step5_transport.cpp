// Testovi za korak 5: transport sedimenta.
//   - seq::advectSediment        (Mei et al., semi-Lagrangian)
//   - seq::transportSedimentFlux (konzervativni, kroz virtualne cijevi; zadani)

#include "TestFramework.h"

#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/Kernels.h"
#include "erosion/Params.h"
#include "erosion/Simulation.h"
#include "erosion/TerrainGen.h"

#include <algorithm>
#include <cfloat>
#include <cmath>

using namespace erosion;

namespace {

double relErr(double actual, double expected)
{
    return std::abs(actual - expected) / std::abs(expected);
}

/// Mreža sa sedimentom s = x + 10*y (linearno polje, bilinearna interpolacija je tačna).
Grid linearSediment(int W, int H, float cellSize)
{
    Grid g(W, H, cellSize);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
            g.s[g.idx(x, y)] = static_cast<float>(x) + 10.0f * static_cast<float>(y);
    return g;
}

/// Stanje nakon nekoliko koraka toka vode preko fBm terena, sa sedimentom.
Grid flowingState(int n)
{
    Grid g(n, n, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 40.0f);
    fillConstant(g.d, 0.5f);
    Params p;
    p.erosionEnabled = false;
    Simulation sim(std::move(g), p);
    sim.run(50); // fluks je sada razvijen
    Grid r = sim.grid();
    for (int y = 0; y < n; ++y)
        for (int x = 0; x < n; ++x)
            r.s[r.idx(x, y)] = 0.1f + 0.05f * static_cast<float>((x * 3 + y * 5) % 7);
    return r;
}

} // namespace

// ── Semi-Lagrangian (Mei) ────────────────────────────────────────────────────

TEST(sl_bez_brzine_sediment_ostaje)
{
    Grid g = linearSediment(16, 12, 2.0f);
    const Field s0 = g.s;
    Field scratch(g.cellCount());
    seq::advectSediment(g, Params{}, scratch);
    CHECK(maxAbsDiff(g.s, s0) == 0.0f);
}

TEST(sl_pomak_za_cijelu_celiju)
{
    // u = l / dt -> pomak tačno jedna ćelija udesno: s_new(x) = s(x - 1).
    Params p;
    Grid g = linearSediment(16, 8, 2.0f);
    fillConstant(g.u, g.cellSize / p.dt);
    const Field s0 = g.s;
    Field scratch(g.cellCount());
    seq::advectSediment(g, p, scratch);

    for (int y = 0; y < 8; ++y) {
        for (int x = 1; x < 16; ++x)
            CHECK_NEAR(g.s[g.idx(x, y)], s0[g.idx(x - 1, y)], 1e-5);
        CHECK_NEAR(g.s[g.idx(0, y)], s0[g.idx(0, y)], 1e-5); // na rubu odsijecanje
    }
}

TEST(sl_pola_celije_daje_prosjek)
{
    // v = -0.5 l / dt -> uzima se iz y + 0.5: prosjek reda y i y + 1.
    Params p;
    Grid g = linearSediment(8, 8, 1.0f);
    fillConstant(g.v, -0.5f * g.cellSize / p.dt);
    const Field s0 = g.s;
    Field scratch(g.cellCount());
    seq::advectSediment(g, p, scratch);
    for (int y = 0; y < 7; ++y)
        for (int x = 0; x < 8; ++x)
            CHECK_NEAR(g.s[g.idx(x, y)], 0.5f * (s0[g.idx(x, y)] + s0[g.idx(x, y + 1)]), 1e-5);
}

TEST(sl_konstantno_polje_ostaje_konstantno)
{
    // Bilinearna interpolacija konstante je ta ista konstanta, za bilo koju brzinu,
    // uključujući ogromne brzine koje izlaze van mreže.
    Grid g(32, 32, 2.0f);
    fillConstant(g.s, 0.7f);
    for (int y = 0; y < 32; ++y)
        for (int x = 0; x < 32; ++x) {
            g.u[g.idx(x, y)] = 500.0f * std::sin(0.4f * x + 0.1f * y);
            g.v[g.idx(x, y)] = 500.0f * std::cos(0.3f * y - 0.2f * x);
        }
    Field scratch(g.cellCount());
    seq::advectSediment(g, Params{}, scratch);
    CHECK(fieldRange(g.s).min >= 0.7f - 1e-6f);
    CHECK(fieldRange(g.s).max <= 0.7f + 1e-6f);
}

// ── Konzervativni transport kroz cijevi ──────────────────────────────────────

TEST(flux_transport_cuva_ukupni_sediment)
{
    Grid g = flowingState(64);
    const double s0 = fieldSum(g.s);
    Field scratch(g.cellCount()), rate(g.cellCount());
    seq::transportSedimentFlux(g, Params{}, scratch, rate);
    CHECK(relErr(fieldSum(g.s), s0) < 1e-6);
    CHECK(fieldRange(g.s).min >= 0.0f);
}

TEST(flux_transport_sediment_prati_vodu)
{
    // Ista koncentracija svuda (s = c * d1): nakon transporta i ažuriranja vode
    // mora vrijediti s = c * d2. Ovo provjerava da sediment ide istim putem kao voda.
    Grid g = flowingState(48);
    Params p;
    seq::computeFlux(g, p); // svjež fluks za trenutno d (= d1)
    const float c = 0.25f;
    for (std::size_t i = 0; i < g.cellCount(); ++i)
        g.s[i] = c * g.d[i];

    Field scratch(g.cellCount()), rate(g.cellCount());
    seq::transportSedimentFlux(g, p, scratch, rate);
    seq::updateWaterAndVelocity(g, p);

    float maxErr = 0.0f;
    for (std::size_t i = 0; i < g.cellCount(); ++i)
        maxErr = std::max(maxErr, std::abs(g.s[i] - c * g.d[i]));
    CHECK(maxErr < 1e-4f);
}

TEST(flux_transport_bez_toka_nista_se_ne_mijenja)
{
    Grid g(16, 16, 2.0f);
    fillConstant(g.d, 1.0f);
    fillConstant(g.s, 0.3f);
    Field scratch(g.cellCount()), rate(g.cellCount());
    seq::transportSedimentFlux(g, Params{}, scratch, rate);
    CHECK(fieldRange(g.s).min == 0.3f && fieldRange(g.s).max == 0.3f);
}

TEST(flux_transport_denormalni_fluks_ne_daje_nan)
{
    // Regresija: fluks reda 1e-40 je ranije dao r = inf i 0 * inf = NaN.
    Grid g(4, 1, 1.0f);
    fillConstant(g.s, 0.5f);
    g.d[1]  = 0.0f;
    g.fR[1] = 1e-40f; // denormalan odlazni fluks iz ćelije bez vode
    g.d[2]  = 1e-30f;
    g.fL[2] = 1e-40f;
    Field scratch(g.cellCount()), rate(g.cellCount());
    seq::transportSedimentFlux(g, Params{}, scratch, rate);
    CHECK(allFinite(g.s));
    CHECK(relErr(fieldSum(g.s), 2.0) < 1e-6);
}

// ── Cijela simulacija ────────────────────────────────────────────────────────

TEST(simulacija_s_transportom_cuva_materijal)
{
    Grid g(128, 128, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);
    const double m0 = totalMaterial(g);

    Params p;
    p.rainRate = 0.02f;
    Simulation sim(std::move(g), p);
    sim.run(3000);

    const Grid& r = sim.grid();
    CHECK(relErr(totalMaterial(r), m0) < 1e-5);
    CHECK(allFinite(r.b) && allFinite(r.s) && allFinite(r.d));
    CHECK(fieldRange(r.s).min >= 0.0f);
}

TEST(sediment_se_talozi_nizvodno_u_bazenu)
{
    // Kosina (x < 64) završava ravnim bazenom (x >= 64). Sediment erodiran na
    // kosini voda nosi dolje, gdje uspori i taloži ga.
    const int W = 96, H = 32;
    Grid g(W, H, 1.0f);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
            g.b[g.idx(x, y)] = 0.3f * static_cast<float>(std::max(0, 64 - x));
    const Field  b0 = g.b;
    const double m0 = totalMaterial(g);

    Simulation sim(std::move(g));
    sim.addSource({4, H / 2, 3.0f, 3.0f});
    sim.run(3000);

    const Grid& r = sim.grid();
    double erodedSlope = 0.0, depositedBasin = 0.0;
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
            const double db = double(r.b[r.idx(x, y)]) - b0[r.idx(x, y)];
            if (x >= 8 && x < 56) erodedSlope   += std::min(0.0, db);
            if (x >= 64)          depositedBasin += std::max(0.0, db);
        }
    CHECK(erodedSlope < -1e-2);
    CHECK(depositedBasin > 1e-3);
    CHECK(relErr(totalMaterial(r), m0) < 1e-5);
}

TEST(semi_lagrangian_gubi_masu_flux_ne)
{
    // Dokumentuje razliku metoda (vidi README, korak 5): ista simulacija,
    // Mei-jeva advekcija gubi materijal, konzervativni transport ne.
    auto run = [](SedimentTransport t) {
        Grid g(96, 96, 2.0f);
        generateFbmTerrain(g, FbmSettings{}, 80.0f);
        const double m0 = totalMaterial(g);
        Params p;
        p.rainRate          = 0.02f;
        p.sedimentTransport = t;
        Simulation sim(std::move(g), p);
        sim.run(1500);
        return (totalMaterial(sim.grid()) - m0) / m0;
    };
    const double lossSL   = run(SedimentTransport::SemiLagrangian);
    const double lossFlux = run(SedimentTransport::Flux);
    CHECK(lossSL < -1e-2);              // gubi više od 1% materijala
    CHECK(std::abs(lossFlux) < 1e-5);   // čuva materijal
}
