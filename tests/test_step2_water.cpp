// Testovi za korak 2: dodavanje vode, fluks kroz virtualne cijevi, visina vode i brzina.

#include "TestFramework.h"

#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/Scenarios.h"
#include "erosion/Simulation.h"
#include "erosion/TerrainGen.h"

#include <algorithm>
#include <cmath>

using namespace erosion;

namespace {

float maxAbs(const Field& f)
{
    float m = 0.0f;
    for (float x : f)
        m = std::max(m, std::abs(x));
    return m;
}

double relErr(double actual, double expected)
{
    return std::abs(actual - expected) / std::abs(expected);
}

Grid fbmGrid(int n, float heightScale, float cellSize = 2.0f)
{
    Grid g(n, n, cellSize);
    generateFbmTerrain(g, FbmSettings{}, heightScale);
    return g;
}

} // namespace

// ── Mirno jezero ─────────────────────────────────────────────────────────────

TEST(mirno_jezero_na_ravnom_terenu_ostaje_mirno)
{
    Grid g(64, 64, 2.0f);
    setupFlatLake(g, 10.0f, 2.0f);
    Simulation sim(std::move(g));
    sim.run(500);

    const Grid& r = sim.grid();
    CHECK(fieldRange(r.d).min == 2.0f);
    CHECK(fieldRange(r.d).max == 2.0f);
    CHECK(maxAbs(r.u) == 0.0f);
    CHECK(maxAbs(r.v) == 0.0f);
    CHECK(maxAbs(r.fL) + maxAbs(r.fR) + maxAbs(r.fT) + maxAbs(r.fB) == 0.0f);
}

TEST(mirno_jezero_na_neravnom_terenu_ostaje_mirno)
{
    // Voda napunjena do ravnog nivoa preko fBm terena: dio ćelija je suh (otoci),
    // ali površina vode je ravna, pa ne smije biti nikakvog toka.
    Grid g = fbmGrid(128, 80.0f);
    fillToLevel(g, 30.0f);
    const Field d0 = g.d;

    Simulation sim(std::move(g));
    sim.run(500);

    CHECK(maxAbsDiff(sim.grid().d, d0) < 1e-5f);
    CHECK(maxAbs(sim.grid().u) < 1e-5f);
    CHECK(maxAbs(sim.grid().v) < 1e-5f);
}

// ── Očuvanje mase i pozitivnost ──────────────────────────────────────────────

TEST(masa_vode_se_cuva_na_fbm_terenu)
{
    Grid g = fbmGrid(128, 80.0f);
    fillConstant(g.d, 1.0f); // sloj vode prati teren, pa voda teče nizbrdo
    const double m0 = totalWater(g);

    Simulation sim(std::move(g));
    float maxCfl = 0.0f;
    for (int i = 0; i < 1000; ++i) {
        sim.step();
        maxCfl = std::max(maxCfl, velocityCfl(sim.grid(), sim.params()).cfl);
    }

    CHECK(relErr(totalWater(sim.grid()), m0) < 1e-4);
    CHECK(fieldRange(sim.grid().d).min >= 0.0f);
    CHECK(allFinite(sim.grid().d));
    CHECK(allFinite(sim.grid().u));
    CHECK(allFinite(sim.grid().v));
    // CFL je po konstrukciji <= 1: K-skaliranje dozvoljava da ćelija u jednom
    // koraku pošalje najviše svu svoju vodu, a tada je u = l / dt (CFL = 1).
    // To se dešava na tankom sloju vode na strmini i nije nestabilnost.
    CHECK(maxCfl <= 1.0f + 1e-4f);
}

TEST(tanak_sloj_vode_na_strmom_terenu_nikad_negativan)
{
    // Najteži slučaj za K-skaliranje: malo vode, velike razlike visina.
    Grid g = fbmGrid(64, 80.0f, 1.0f);
    fillConstant(g.d, 0.01f);
    const double m0 = totalWater(g);

    Simulation sim(std::move(g));
    for (int i = 0; i < 500; ++i) {
        sim.step();
        if (fieldRange(sim.grid().d).min < 0.0f) {
            CHECK(false);
            break;
        }
    }
    CHECK(relErr(totalWater(sim.grid()), m0) < 1e-4);
}

TEST(rubovi_su_zatvoreni)
{
    Grid g = fbmGrid(48, 40.0f);
    fillConstant(g.d, 2.0f);
    Simulation sim(std::move(g));
    sim.run(100);

    const Grid& r = sim.grid();
    for (int y = 0; y < r.height; ++y) {
        CHECK(r.fL[r.idx(0, y)] == 0.0f);
        CHECK(r.fR[r.idx(r.width - 1, y)] == 0.0f);
    }
    for (int x = 0; x < r.width; ++x) {
        CHECK(r.fT[r.idx(x, 0)] == 0.0f);
        CHECK(r.fB[r.idx(x, r.height - 1)] == 0.0f);
    }
}

// ── Dodavanje vode ───────────────────────────────────────────────────────────

TEST(kisa_dodaje_tacno_ocekivanu_kolicinu)
{
    Grid g = fbmGrid(64, 40.0f);
    Params p;
    p.rainRate = 0.5f;
    Simulation sim(std::move(g), p);
    sim.run(100);

    const double area     = 2.0 * 2.0;
    const double expected = 0.5 * p.dt * 100 * 64 * 64 * area;
    CHECK(relErr(totalWater(sim.grid()), expected) < 1e-5);
}

TEST(izvor_vode_dodaje_samo_unutar_kruga)
{
    // Samo faza dodavanja vode (bez toka), da voda ne stigne oteći iz kruga.
    Grid              g(41, 41, 1.0f);
    const Params      p;
    const WaterSource src{20, 20, 3.0f, 2.0f};
    seq::addWater(g, p, {&src, 1});

    // Tačke (dx, dy) s dx² + dy² <= 9: ima ih 29.
    const double expected = 29 * 2.0 * p.dt;
    CHECK(relErr(totalWater(g), expected) < 1e-5);
    CHECK(g.d[g.idx(20, 23)] > 0.0f);  // na rubu kruga (udaljenost 3)
    CHECK(g.d[g.idx(20, 24)] == 0.0f); // izvan kruga (udaljenost 4)
    CHECK(g.d[g.idx(23, 22)] == 0.0f); // izvan kruga (9 + 4 > 9)
}

// ── Dam-break ────────────────────────────────────────────────────────────────

TEST(dam_break_je_simetrican_i_cuva_masu)
{
    const int N = 65; // neparno, da centar padne tačno na ćeliju
    Grid g(N, N, 1.0f);
    setupDamBreak(g, 1.0f, 3.0f, 8);
    const double m0 = totalWater(g);

    Simulation sim(std::move(g));
    const int c        = N / 2;
    float     ringPeak = 0.0f; // najveća dubina u tački 20 ćelija od centra
    float     maxCfl   = 0.0f;
    for (int i = 0; i < 300; ++i) {
        sim.step();
        ringPeak = std::max(ringPeak, sim.grid().d[sim.grid().idx(c + 20, c)]);
        maxCfl   = std::max(maxCfl, velocityCfl(sim.grid(), sim.params()).cfl);
    }

    const Grid& r = sim.grid();
    CHECK(relErr(totalWater(r), m0) < 1e-5);
    CHECK(fieldRange(r.d).min >= 0.0f);
    CHECK(maxCfl < 1.0f);

    // Stub se srušio, a talas je stigao do prstena (podigao vodu iznad 1).
    CHECK(r.d[r.idx(c, c)] < 3.0f);
    CHECK(ringPeak > 1.05f);

    // Simetrija: ogledalo po x, po y i po dijagonali (zamjena x i y).
    float asymX = 0.0f, asymY = 0.0f, asymT = 0.0f;
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            const float h = r.d[r.idx(x, y)];
            asymX = std::max(asymX, std::abs(h - r.d[r.idx(N - 1 - x, y)]));
            asymY = std::max(asymY, std::abs(h - r.d[r.idx(x, N - 1 - y)]));
            asymT = std::max(asymT, std::abs(h - r.d[r.idx(y, x)]));
        }
    }
    CHECK(asymX < 1e-4f);
    CHECK(asymY < 1e-4f);
    CHECK(asymT < 1e-4f);
}

// ── Smjer toka i brzina ──────────────────────────────────────────────────────

TEST(voda_tece_nizbrdo)
{
    // Teren raste s x, pa voda mora teći prema x = 0: u < 0.
    const int W = 64, H = 16;
    Grid g(W, H, 1.0f);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
            g.b[g.idx(x, y)] = 0.2f * static_cast<float>(x);
    fillConstant(g.d, 0.5f);

    auto centroidX = [](const Grid& gr) {
        double sum = 0.0, wsum = 0.0;
        for (int y = 0; y < gr.height; ++y)
            for (int x = 0; x < gr.width; ++x) {
                sum  += x * static_cast<double>(gr.d[gr.idx(x, y)]);
                wsum += gr.d[gr.idx(x, y)];
            }
        return sum / wsum;
    };
    const double cx0 = centroidX(g);

    Simulation sim(std::move(g));
    sim.run(20);
    // Na sredini kosine (dalje od zida) brzina mora biti negativna.
    CHECK(sim.grid().u[sim.grid().idx(W / 2, H / 2)] < 0.0f);
    CHECK(std::abs(sim.grid().v[sim.grid().idx(W / 2, H / 2)]) < 1e-5f); // nema toka po y

    sim.run(2000);
    CHECK(centroidX(sim.grid()) < cx0 - 10.0); // voda se skupila na nižoj strani
}

TEST(simulacija_broji_korake_i_vrijeme)
{
    Simulation sim(Grid(8, 8, 1.0f));
    sim.run(10);
    CHECK(sim.stepCount() == 10);
    CHECK_NEAR(sim.time(), 10 * sim.params().dt, 1e-6);
    CHECK(sim.timings().steps == 10);
    sim.resetTimings();
    CHECK(sim.timings().steps == 0);
}
