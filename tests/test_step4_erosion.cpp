// Testovi za korak 4: erozija i depozicija.

#include "TestFramework.h"

#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/Kernels.h"
#include "erosion/Params.h"
#include "erosion/Scenarios.h"
#include "erosion/Simulation.h"
#include "erosion/TerrainGen.h"

#include <algorithm>
#include <cmath>

using namespace erosion;

namespace {

double relErr(double actual, double expected)
{
    return std::abs(actual - expected) / std::abs(expected);
}

/// 3 x 3 mreža, l = 1, teren b = 2x (nagib 2 po x, 0 po y).
/// Centralna ćelija: sin α = 2 / sqrt(5), brzina (3, 4) -> |v| = 5, d = 1.
Grid slopeCell(float s0, float depth = 1.0f)
{
    Grid g(3, 3, 1.0f);
    for (int y = 0; y < 3; ++y)
        for (int x = 0; x < 3; ++x)
            g.b[g.idx(x, y)] = 2.0f * static_cast<float>(x);
    const std::size_t c = g.idx(1, 1);
    g.u[c] = 3.0f;
    g.v[c] = 4.0f;
    g.d[c] = depth;
    g.s[c] = s0;
    return g;
}

const double kSinA = 2.0 / std::sqrt(5.0);

} // namespace

// ── Formula na jednoj ćeliji ─────────────────────────────────────────────────

TEST(erozija_kad_je_kapacitet_veci_od_sedimenta)
{
    Params p; // Kc = 0.5, Ks = 0.5/s, dt = 0.02 -> ks = 0.01
    Grid   g = slopeCell(0.5f);
    Field  scratch(g.cellCount());
    seq::erodeDeposit(g, p, scratch);

    const double C      = 0.5 * kSinA * 5.0;     // 2.236068
    const double amount = 0.01 * (C - 0.5);      // 0.01736068
    CHECK_NEAR(g.b[g.idx(1, 1)], 2.0 - amount, 1e-6);
    CHECK_NEAR(g.s[g.idx(1, 1)], 0.5 + amount, 1e-6);
}

TEST(depozicija_kad_je_sediment_veci_od_kapaciteta)
{
    Params p; // Kd = 1/s, dt = 0.02 -> kd = 0.02
    Grid   g = slopeCell(3.0f);
    Field  scratch(g.cellCount());
    seq::erodeDeposit(g, p, scratch);

    const double C      = 0.5 * kSinA * 5.0;
    const double amount = 0.02 * (3.0 - C);
    CHECK_NEAR(g.b[g.idx(1, 1)], 2.0 + amount, 1e-6);
    CHECK_NEAR(g.s[g.idx(1, 1)], 3.0 - amount, 1e-6);
}

TEST(plitka_voda_smanjuje_kapacitet_lmax)
{
    Params p;
    p.maxErosionDepth = 0.1f;
    Grid   g = slopeCell(0.0f, 0.05f); // d = Kdmax / 2 -> lmax = 0.5
    Field  scratch(g.cellCount());
    seq::erodeDeposit(g, p, scratch);

    const double C = 0.5 * kSinA * 5.0 * 0.5;
    CHECK_NEAR(g.s[g.idx(1, 1)], 0.01 * C, 1e-6);
}

TEST(suha_celija_se_ne_erodira)
{
    // Ćelije oko centra imaju d = 0, brzinu 0 i s = 0: ne smiju se promijeniti.
    Params p;
    Grid   g = slopeCell(0.0f);
    Field  scratch(g.cellCount());
    seq::erodeDeposit(g, p, scratch);
    CHECK(g.b[g.idx(0, 0)] == 0.0f);
    CHECK(g.b[g.idx(2, 2)] == 4.0f);
    CHECK(g.s[g.idx(2, 1)] == 0.0f);
}

// ── Cijela mreža naspram referentne implementacije ───────────────────────────

TEST(erozija_odgovara_referentnoj_implementaciji)
{
    // Proizvoljno, ali determinističko stanje: fBm teren, valovita voda,
    // brzina i sediment. Referenca računa u double i uvijek čita STARI teren.
    const int N = 48;
    Grid g(N, N, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 40.0f);
    for (int y = 0; y < N; ++y)
        for (int x = 0; x < N; ++x) {
            const std::size_t i = g.idx(x, y);
            g.d[i] = std::max(0.0f, 0.3f * std::sin(0.3f * x) + 0.2f * std::cos(0.2f * y));
            g.u[i] = 2.0f * std::sin(0.17f * (x + y));
            g.v[i] = 1.5f * std::cos(0.23f * (x - y));
            g.s[i] = 0.05f * static_cast<float>((x * 7 + y * 13) % 11);
        }

    const Grid before = g;
    Params     p;
    Field      scratch(g.cellCount());
    seq::erodeDeposit(g, p, scratch);

    double maxErrB = 0.0, maxErrS = 0.0;
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            const int xl = std::max(x - 1, 0), xr = std::min(x + 1, N - 1);
            const int yu = std::max(y - 1, 0), yd = std::min(y + 1, N - 1);
            const double gx = (double(before.b[before.idx(xr, y)]) - before.b[before.idx(xl, y)]) / ((xr - xl) * 2.0);
            const double gy = (double(before.b[before.idx(x, yd)]) - before.b[before.idx(x, yu)]) / ((yd - yu) * 2.0);
            const double g2   = gx * gx + gy * gy;
            const double sinA = std::max<double>(p.minTilt, std::sqrt(g2 / (1.0 + g2)));

            const std::size_t i = before.idx(x, y);
            const double speed = std::hypot(double(before.u[i]), double(before.v[i]));
            const double lmax  = std::min(1.0, before.d[i] / double(p.maxErosionDepth));
            const double C     = p.sedimentCapacity * sinA * speed * lmax;
            const double s0    = before.s[i];
            const double amt   = C > s0 ? double(p.dissolvingRate) * p.dt * (C - s0)
                                        : -double(p.depositionRate) * p.dt * (s0 - C);
            maxErrB = std::max(maxErrB, std::abs(g.b[i] - (before.b[i] - amt)));
            maxErrS = std::max(maxErrS, std::abs(g.s[i] - (s0 + amt)));
        }
    }
    CHECK(maxErrB < 1e-5);
    CHECK(maxErrS < 1e-5);

    // Nakon zamjene bafera, scratch sadrži stari teren.
    CHECK(maxAbsDiff(scratch, before.b) == 0.0f);
}

// ── Cijela simulacija ────────────────────────────────────────────────────────

TEST(ukupni_materijal_se_cuva)
{
    Grid g(128, 128, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);
    fillConstant(g.d, 1.0f);
    const double m0 = totalMaterial(g);

    Params p;
    p.rainRate = 0.02f;
    Simulation sim(std::move(g), p);
    sim.run(2000);

    const Grid& r = sim.grid();
    CHECK(relErr(totalMaterial(r), m0) < 1e-5);
    CHECK(fieldRange(r.s).min >= 0.0f);
    CHECK(allFinite(r.b));
    CHECK(allFinite(r.s));
    // Erozija se stvarno desila: teren se promijenio.
    CHECK(fieldRange(r.s).max > 0.0f);
}

TEST(mirno_jezero_ne_erodira)
{
    Grid g(96, 96, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);
    fillToLevel(g, 30.0f);
    const Field b0 = g.b;

    // Bez isparavanja: ono bi spustilo nivo različito u plitkim i dubokim
    // ćelijama, pa površina više ne bi bila ravna i voda bi potekla.
    Params p;
    p.evaporationRate = 0.0f;
    Simulation sim(std::move(g), p);
    sim.run(500);
    CHECK(maxAbsDiff(sim.grid().b, b0) == 0.0f);
    CHECK(fieldRange(sim.grid().s).max == 0.0f);
}

TEST(izvor_na_kosini_erodira_kosinu)
{
    // Kosina pada prema x = W-1 i završava ravnim bazenom.
    // Napomena: taloženje NIZVODNO (u bazenu) testiramo tek u koraku 5, jer
    // bez advekcije sediment ostaje u ćeliji u kojoj je erodiran.
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
    double erodedSlope = 0.0, sedimentSlope = 0.0;
    for (int y = 0; y < H; ++y)
        for (int x = 8; x < 56; ++x) {
            erodedSlope   += std::min(0.0, double(r.b[r.idx(x, y)]) - b0[r.idx(x, y)]);
            sedimentSlope += r.s[r.idx(x, y)];
        }
    CHECK(erodedSlope < -1e-2);  // kosina je izgubila materijal
    CHECK(sedimentSlope > 1e-2); // i on je sada u vodi kao sediment
    CHECK(relErr(totalMaterial(r), m0) < 1e-5);
}

TEST(erozija_se_moze_iskljuciti)
{
    Grid g(64, 64, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);
    fillConstant(g.d, 1.0f);
    const Field b0 = g.b;

    Params p;
    p.erosionEnabled = false;
    Simulation sim(std::move(g), p);
    sim.run(200);
    CHECK(maxAbsDiff(sim.grid().b, b0) == 0.0f);
    CHECK(fieldRange(sim.grid().s).max == 0.0f);
}
