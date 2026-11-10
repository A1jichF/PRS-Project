// Testovi za stijenu (bedrock): polje Grid::r i sloj tla iznad njega.

#include "TestFramework.h"

#include "erosion/Bedrock.h"
#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/Kernels.h"
#include "erosion/Params.h"
#include "erosion/Simulation.h"
#include "erosion/TerrainGen.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>

using namespace erosion;

namespace {

std::uintptr_t pageOffset(const Field& f)
{
    return reinterpret_cast<std::uintptr_t>(f.data()) % 4096;
}

/// true ako nijedna dva polja mreže nisu na istom pomaku unutar 4 KB stranice.
bool allOffsetsDistinct(const Grid& g)
{
    const Field* fields[] = {&g.b, &g.d, &g.s, &g.fL, &g.fR, &g.fT, &g.fB, &g.u, &g.v, &g.r};
    for (std::size_t a = 0; a < std::size(fields); ++a)
        for (std::size_t b = a + 1; b < std::size(fields); ++b)
            if (pageOffset(*fields[a]) == pageOffset(*fields[b]))
                return false;
    return true;
}

} // namespace

TEST(stijena_clear_vraca_na_bez_stijene)
{
    Grid g(16, 8, 1.0f);
    fillConstant(g.r, 3.0f);
    fillConstant(g.b, 5.0f);
    g.clear();
    CHECK(fieldRange(g.r).min == kNoBedrock);
    CHECK(fieldRange(g.r).max == kNoBedrock);
    CHECK(fieldSum(g.b) == 0.0);
}

TEST(stijena_resize_vraca_na_bez_stijene)
{
    Grid g(16, 8, 1.0f);
    fillConstant(g.r, 3.0f);
    g.resize(20, 10, 2.0f);
    CHECK(g.r.size() == g.cellCount());
    CHECK(fieldRange(g.r).max == kNoBedrock);
}

TEST(stijena_bez_stijene_je_konacna_vrijednost)
{
    // kNoBedrock nije -inf, pa provjere konačnosti i dalje rade nad cijelom mrežom.
    Grid g(8, 8, 1.0f);
    CHECK(allFinite(g.r));
    // Debljina tla b - r mora biti ogromna i konačna za svaku realnu visinu terena.
    CHECK(1000.0f - kNoBedrock > 1e30f);
    CHECK(-1000.0f - kNoBedrock > 1e30f);
}

TEST(stijena_se_ne_poklapa_s_pomocnim_baferima_simulacije)
{
    // Zaštita od 4K aliasinga: pomoćni baferi simulacije moraju biti na
    // pomacima različitim od SVIH polja mreže, uključujući r (koje se u
    // kernelu erozije koristi zajedno s pomoćnim baferom terena).
    //
    // Pomoćni baferi su privatni, ali kerneli ih zamjenjuju s poljima mreže:
    // poslije jednog koraka g.b je bivši bafer terena, a g.s bivši bafer
    // sedimenta. Poslije dva koraka su opet originalni.
    Grid g(64, 64, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);
    fillConstant(g.d, 0.5f);

    Simulation sim(std::move(g), Params{});
    CHECK(allOffsetsDistinct(sim.grid()));

    const std::uintptr_t bBefore = pageOffset(sim.grid().b);
    const std::uintptr_t sBefore = pageOffset(sim.grid().s);
    sim.step();
    // Zamjena se zaista desila, pa provjera ispod gleda pomoćne bafere.
    CHECK(pageOffset(sim.grid().b) != bBefore);
    CHECK(pageOffset(sim.grid().s) != sBefore);
    CHECK(allOffsetsDistinct(sim.grid()));

    sim.step();
    CHECK(pageOffset(sim.grid().b) == bBefore);
    CHECK(pageOffset(sim.grid().s) == sBefore);
}

TEST(stijena_simulacija_je_ne_mijenja)
{
    // Uz zadanu tvrdoću (rockErodibility = 0) simulacija nikad ne mijenja r.
    Grid g(48, 48, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);
    for (std::size_t i = 0; i < g.cellCount(); ++i)
        g.r[i] = g.b[i] - 1.5f;
    const Field r0 = g.r;

    Params p;
    p.rainRate = 0.02f;
    Simulation sim(std::move(g), p);
    sim.run(100);
    CHECK(maxAbsDiff(sim.grid().r, r0) == 0.0f);
}

// ── Kernel erozije sa stijenom: jedna ćelija ─────────────────────────────────
//
// Ista ćelija kao u test_step4_erosion.cpp: 3 x 3, l = 1, teren b = 2x, u
// centru brzina (3, 4), d = 1. Kapacitet C = 0.5 * (2/sqrt(5)) * 5 = 2.236068,
// ks = Ks*dt = 0.01, pa uz s = 0.5 erozija "želi" 0.01 * (C - 0.5) = 0.01736068.

namespace {

constexpr double kWanted = 0.01736068;

Grid slopeCell(float s0)
{
    Grid g(3, 3, 1.0f);
    for (int y = 0; y < 3; ++y)
        for (int x = 0; x < 3; ++x)
            g.b[g.idx(x, y)] = 2.0f * static_cast<float>(x);
    const std::size_t c = g.idx(1, 1);
    g.u[c] = 3.0f;
    g.v[c] = 4.0f;
    g.d[c] = 1.0f;
    g.s[c] = s0;
    return g;
}

/// Jedan poziv kernela erozije; vraća mrežu poslije.
Grid erodeOnce(Grid g, const Params& p)
{
    Field scratch(g.cellCount());
    seq::erodeDeposit(g, p, scratch);
    return g;
}

double relErr(double actual, double expected)
{
    return std::abs(actual - expected) / std::abs(expected);
}

} // namespace

TEST(stijena_gola_stijena_se_ne_erodira)
{
    Grid g = slopeCell(0.5f);
    const std::size_t c = g.idx(1, 1);
    g.r[c] = g.b[c]; // nema tla

    const Grid out = erodeOnce(g, Params{});
    CHECK(out.b[c] == 2.0f);
    CHECK(out.r[c] == 2.0f);
    CHECK(out.s[c] == 0.5f);
}

TEST(stijena_tanko_tlo_se_odnese_tacno_do_stijene)
{
    Grid g = slopeCell(0.5f);
    const std::size_t c = g.idx(1, 1);
    g.r[c] = 1.995f; // 0.005 tla, a erozija želi 0.0174

    const Grid out = erodeOnce(g, Params{});
    CHECK(out.b[c] == 1.995f); // tačno na stijeni, ne ispod
    CHECK(out.r[c] == 1.995f);
    CHECK_NEAR(out.s[c], 0.5 + 0.005, 1e-6);
}

TEST(stijena_duboko_tlo_ne_mijenja_eroziju)
{
    // Dok ima dovoljno tla, rezultat je bit-identičan onome bez stijene.
    Grid g = slopeCell(0.5f);
    const std::size_t c = g.idx(1, 1);
    const Grid noRock = erodeOnce(g, Params{});

    g.r[c] = 1.0f; // 1 m tla
    const Grid out = erodeOnce(g, Params{});
    CHECK(out.b[c] == noRock.b[c]);
    CHECK(out.s[c] == noRock.s[c]);
    CHECK(out.r[c] == 1.0f);
    CHECK_NEAR(out.b[c], 2.0 - kWanted, 1e-6);
}

TEST(stijena_tvrdoca_usporava_eroziju_stijene)
{
    Grid g = slopeCell(0.5f);
    const std::size_t c = g.idx(1, 1);
    g.r[c] = g.b[c];

    Params p;
    p.rockErodibility = 0.25f;
    const Grid out = erodeOnce(g, p);
    CHECK_NEAR(out.r[c], 2.0 - 0.25 * kWanted, 1e-6);
    CHECK(out.b[c] == out.r[c]); // stijena ostaje gola
    CHECK_NEAR(out.s[c], 0.5 + 0.25 * kWanted, 1e-6);
}

TEST(stijena_prvo_tlo_pa_stijena)
{
    // 0.005 tla, tvrdoća 0.5: cijelo tlo + polovina ostatka iz stijene.
    Grid g = slopeCell(0.5f);
    const std::size_t c = g.idx(1, 1);
    g.r[c] = 1.995f;

    Params p;
    p.rockErodibility = 0.5f;
    const Grid   out      = erodeOnce(g, p);
    const double fromRock = 0.5 * (kWanted - 0.005);
    CHECK_NEAR(out.r[c], 1.995 - fromRock, 1e-6);
    CHECK(out.b[c] == out.r[c]);
    CHECK_NEAR(out.s[c], 0.5 + 0.005 + fromRock, 1e-6);
}

TEST(stijena_tvrdoca_se_ogranicava_na_0_1)
{
    Grid g = slopeCell(0.5f);
    const std::size_t c = g.idx(1, 1);
    g.r[c] = g.b[c];

    Params p;
    p.rockErodibility = 1.0f;
    const Grid asSoil = erodeOnce(g, p);
    p.rockErodibility = 7.0f;
    CHECK(erodeOnce(g, p).b[c] == asSoil.b[c]);
    p.rockErodibility = -3.0f;
    CHECK(erodeOnce(g, p).b[c] == 2.0f);

    // Tvrdoća 1: stijena se erodira jednako kao tlo.
    CHECK_NEAR(asSoil.b[c], 2.0 - kWanted, 1e-6);
}

TEST(stijena_talog_na_stijeni_se_opet_moze_odnijeti)
{
    // s > C: sediment se taloži na golu stijenu kao novo tlo, r se ne mijenja.
    Grid g = slopeCell(3.0f);
    const std::size_t c = g.idx(1, 1);
    g.r[c] = g.b[c];

    Field scratch(g.cellCount());
    seq::erodeDeposit(g, Params{}, scratch);
    CHECK(g.b[c] > 2.0f);
    CHECK(g.r[c] == 2.0f);

    // Čista voda zatim odnese talog, ali samo do stijene.
    for (int i = 0; i < 5; ++i) {
        g.s[c] = 0.0f;
        seq::erodeDeposit(g, Params{}, scratch);
    }
    CHECK(g.b[c] == 2.0f);
    CHECK(g.r[c] == 2.0f);
}

TEST(stijena_iznad_povrsine_ne_stvara_materijal)
{
    // Neispravan ulaz (r > b): stijena se računa od površine, teren se ne diže.
    Grid g = slopeCell(0.5f);
    const std::size_t c = g.idx(1, 1);
    g.r[c] = 5.0f;

    const Grid out = erodeOnce(g, Params{});
    CHECK(out.b[c] == 2.0f);
    CHECK(out.r[c] == 2.0f);
    CHECK(out.s[c] == 0.5f);
}

// ── Cijela simulacija ────────────────────────────────────────────────────────

namespace {

/// fBm 96 x 96 s 1 m vode i kišom; isti scenarij sa i bez stijene.
Grid rainyTerrain()
{
    Grid g(96, 96, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);
    fillConstant(g.d, 1.0f);
    return g;
}

constexpr float kSoil  = 0.05f;
constexpr int   kSteps = 1500;

} // namespace

TEST(stijena_teren_nikad_ne_ide_ispod_stijene)
{
    Grid g = rainyTerrain();
    const Field b0 = g.b;
    for (std::size_t i = 0; i < g.cellCount(); ++i)
        g.r[i] = g.b[i] - kSoil;
    const Field  r0 = g.r;
    const double m0 = totalMaterial(g);

    Params p;
    p.rainRate = 0.02f;
    Simulation sim(g, p);
    sim.run(kSteps);
    const Grid& out = sim.grid();

    std::size_t bare      = 0;
    bool        aboveRock = true;
    for (std::size_t i = 0; i < out.cellCount(); ++i) {
        aboveRock = aboveRock && out.b[i] >= out.r[i];
        bare += out.b[i] == out.r[i] ? 1 : 0;
    }
    CHECK(aboveRock);
    CHECK(maxAbsDiff(out.r, r0) == 0.0f); // tvrdoća 0: stijena netaknuta
    CHECK(bare > 0);                      // erozija je negdje stvarno došla do stijene
    CHECK(relErr(totalMaterial(out), m0) < 1e-5);

    // Test ima smisla samo ako bi bez stijene erozija išla dublje od sloja tla.
    Simulation noRock(rainyTerrain(), p);
    noRock.run(kSteps);
    float deepest = 0.0f;
    for (std::size_t i = 0; i < b0.size(); ++i)
        deepest = std::min(deepest, noRock.grid().b[i] - b0[i]);
    CHECK(deepest < -2.0f * kSoil);
}

TEST(stijena_meka_stijena_se_spusta_ali_ostaje_ispod_terena)
{
    Grid g = rainyTerrain();
    for (std::size_t i = 0; i < g.cellCount(); ++i)
        g.r[i] = g.b[i] - kSoil;
    const Field  r0 = g.r;
    const double m0 = totalMaterial(g);

    Params p;
    p.rainRate        = 0.02f;
    p.rockErodibility = 0.1f;
    Simulation sim(std::move(g), p);
    sim.run(kSteps);
    const Grid& out = sim.grid();

    bool  aboveRock = true, rockOnlyDown = true;
    float maxRockDrop = 0.0f;
    for (std::size_t i = 0; i < out.cellCount(); ++i) {
        aboveRock    = aboveRock && out.b[i] >= out.r[i];
        rockOnlyDown = rockOnlyDown && out.r[i] <= r0[i];
        maxRockDrop  = std::max(maxRockDrop, r0[i] - out.r[i]);
    }
    CHECK(aboveRock);
    CHECK(rockOnlyDown);
    CHECK(maxRockDrop > 0.0f);
    CHECK(relErr(totalMaterial(out), m0) < 1e-5);
}

TEST(stijena_korito_od_stijene_se_ne_pomjera_dok_voda_tece)
{
    // Kosina od gole stijene s izvorom na vrhu: voda teče, a teren ostaje
    // bit-identičan (nema ni erozije ni sedimenta koji bi se taložio).
    const int W = 96, H = 32;
    Grid g(W, H, 1.0f);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
            g.b[g.idx(x, y)] = 0.3f * static_cast<float>(std::max(0, 64 - x));
    g.r = g.b;
    const Field b0 = g.b;

    Simulation sim(std::move(g));
    sim.addSource({4, H / 2, 3.0f, 3.0f});
    sim.run(3000);
    const Grid& out = sim.grid();

    CHECK(maxAbsDiff(out.b, b0) == 0.0f);
    CHECK(maxAbsDiff(out.r, b0) == 0.0f);
    CHECK(fieldRange(out.s).max == 0.0f);
    CHECK(fieldRange(out.u).max > 0.5f); // voda zaista teče niz kosinu
}

// ── Zadavanje stijene (Bedrock.h) ────────────────────────────────────────────

namespace {

const float kNaN = std::numeric_limits<float>::quiet_NaN();
const float kInf = std::numeric_limits<float>::infinity();

/// 4 x 1 mreža s visinama terena 1, 2, 3, 4 (l = 2).
Grid stairs()
{
    Grid g(4, 1, 2.0f);
    for (int x = 0; x < 4; ++x)
        g.b[g.idx(x, 0)] = static_cast<float>(x + 1);
    return g;
}

} // namespace

TEST(stijena_isti_sloj_tla_svuda)
{
    Grid g = stairs();
    CHECK(!hasBedrock(g));

    setBedrockUniformSoil(g, 0.5f);
    CHECK(hasBedrock(g));
    for (int x = 0; x < 4; ++x)
        CHECK(g.r[g.idx(x, 0)] == g.b[g.idx(x, 0)] - 0.5f);

    setBedrockUniformSoil(g, -3.0f); // negativno = 0 = gola stijena
    CHECK(maxAbsDiff(g.r, g.b) == 0.0f);

    setBedrockUniformSoil(g, kInf); // beskonačno tlo = nema stijene
    CHECK(!hasBedrock(g));
    setBedrockUniformSoil(g, 0.5f);
    setBedrockUniformSoil(g, kNaN);
    CHECK(!hasBedrock(g));
}

TEST(stijena_sloj_tla_po_celiji)
{
    Grid        g       = stairs();
    const float depth[] = {0.0f, 0.25f, -1.0f, kNaN};
    setBedrockFromSoilDepth(g, depth);
    CHECK(g.r[0] == 1.0f);       // gola stijena
    CHECK(g.r[1] == 1.75f);
    CHECK(g.r[2] == 3.0f);       // negativno = gola stijena
    CHECK(g.r[3] == kNoBedrock); // NaN = nema stijene
    CHECK(allFinite(g.r));
}

TEST(stijena_do_zadane_visine)
{
    Grid g = stairs();
    setBedrockLevel(g, 2.5f);
    CHECK(g.r[0] == 1.0f); // ispod nivoa: gola stijena
    CHECK(g.r[1] == 2.0f);
    CHECK(g.r[2] == 2.5f); // iznad nivoa: tlo 0.5
    CHECK(g.r[3] == 2.5f); // tlo 1.5

    setBedrockLevel(g, kNaN);
    CHECK(!hasBedrock(g));
    setBedrockLevel(g, 2.5f);
    setBedrockLevel(g, -kInf);
    CHECK(!hasBedrock(g));
}

TEST(stijena_proizvoljna_mapa_se_ogranicava_na_teren)
{
    Grid        g         = stairs();
    const float heights[] = {0.5f, 9.0f, kInf, 3.5f};
    setBedrockHeights(g, heights);
    CHECK(g.r[0] == 0.5f);
    CHECK(g.r[1] == 2.0f);       // stijena iznad terena se spušta na teren
    CHECK(g.r[2] == kNoBedrock); // beskonačno = nema stijene
    CHECK(g.r[3] == 3.5f);
    for (std::size_t i = 0; i < g.cellCount(); ++i)
        CHECK(g.r[i] <= g.b[i]);

    clearBedrock(g);
    CHECK(!hasBedrock(g));
}

TEST(stijena_statistika_tla)
{
    Grid g = stairs();
    CHECK(bedrockStats(g).rockCells == 0);
    CHECK(bedrockStats(g).soilVolume == 0.0);

    const float depth[] = {0.0f, 0.25f, 1.0f, kNaN};
    setBedrockFromSoilDepth(g, depth);
    const BedrockStats st = bedrockStats(g);
    CHECK(st.rockCells == 3); // ćelija bez stijene se ne broji
    CHECK(st.bareCells == 1);
    CHECK_NEAR(st.soilVolume, (0.25 + 1.0) * 4.0, 1e-9); // l² = 4
    CHECK(bedrockStats(g, 0.3f).bareCells == 2);         // prag za "golo"
}

TEST(stijena_zadana_slojem_tla_ogranicava_simulaciju)
{
    // Isti scenarij kao gore, ali stijena zadana funkcijom iz Bedrock.h:
    // tla ukupno može nestati najviše onoliko koliko ga je bilo.
    Grid g = rainyTerrain();
    setBedrockUniformSoil(g, kSoil);
    const BedrockStats before = bedrockStats(g);
    CHECK(before.bareCells == 0);

    Params p;
    p.rainRate = 0.02f;
    Simulation sim(std::move(g), p);
    sim.run(kSteps);

    const BedrockStats after = bedrockStats(sim.grid());
    CHECK(after.rockCells == before.rockCells);
    CHECK(after.bareCells > 0);
    // Tlo se seli (erozija + taloženje), a dio je u vodi kao sediment.
    const double l2        = 4.0;
    const double suspended = fieldSum(sim.grid().s) * l2;
    CHECK(relErr(after.soilVolume + suspended, before.soilVolume) < 1e-3);
}
