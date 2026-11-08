// Testovi za korak 1: Grid, generator terena, PGM ulaz/izlaz, dijagnostika.

#include "TestFramework.h"

#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/ImageIO.h"
#include "erosion/TerrainGen.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>

using namespace erosion;

namespace {

std::string tempPath(const char* name)
{
    return (std::filesystem::temp_directory_path() / name).string();
}

} // namespace

TEST(grid_alocira_sva_polja_na_nulu)
{
    Grid g(37, 21, 2.0f);
    CHECK(g.cellCount() == 37u * 21u);
    for (const Field* f : {&g.b, &g.d, &g.s, &g.fL, &g.fR, &g.fT, &g.fB, &g.u, &g.v}) {
        CHECK(f->size() == g.cellCount());
        CHECK(fieldSum(*f) == 0.0);
    }
}

TEST(grid_polja_su_poravnata_na_64_bajta)
{
    Grid g(100, 3, 1.0f);
    for (const Field* f : {&g.b, &g.d, &g.s, &g.fL, &g.fR, &g.fT, &g.fB, &g.u, &g.v})
        CHECK(reinterpret_cast<std::uintptr_t>(f->data()) % 64 == 0);
}

TEST(grid_indeks_je_red_po_red)
{
    Grid g(5, 4, 1.0f);
    CHECK(g.idx(0, 0) == 0);
    CHECK(g.idx(4, 0) == 4);
    CHECK(g.idx(0, 1) == 5);  // novi red počinje nakon width ćelija
    CHECK(g.idx(4, 3) == 19); // zadnja ćelija
    CHECK(g.inBounds(4, 3));
    CHECK(!g.inBounds(5, 0));
    CHECK(!g.inBounds(0, -1));
}

TEST(grid_odbija_neispravne_dimenzije)
{
    bool threw = false;
    try { Grid g(0, 10, 1.0f); } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
}

TEST(fbm_je_deterministican)
{
    Grid a(128, 128, 1.0f), b(128, 128, 1.0f);
    FbmSettings s;
    s.seed = 7;
    generateFbmTerrain(a, s, 80.0f);
    generateFbmTerrain(b, s, 80.0f);
    CHECK(maxAbsDiff(a.b, b.b) == 0.0f);

    Grid c(128, 128, 1.0f);
    s.seed = 8;
    generateFbmTerrain(c, s, 80.0f);
    CHECK(maxAbsDiff(a.b, c.b) > 1.0f); // drugi seed daje drugi teren
}

TEST(fbm_teren_je_u_opsegu_0_heightScale)
{
    Grid g(200, 150, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);
    const FieldRange r = fieldRange(g.b);
    CHECK_NEAR(r.min, 0.0f, 1e-4);
    CHECK_NEAR(r.max, 80.0f, 1e-3);
    CHECK(allFinite(g.b));
}

TEST(fbm_ne_dira_ostala_polja)
{
    Grid g(64, 64, 1.0f);
    generateFbmTerrain(g, FbmSettings{}, 10.0f);
    CHECK(fieldSum(g.d) == 0.0);
    CHECK(fieldSum(g.s) == 0.0);
}

TEST(fbm_isti_oblik_na_razlicitim_rezolucijama)
{
    // Ista tačka u normalizovanim koordinatama mora imati sličnu visinu na 129² i 257².
    Grid lo(129, 129, 1.0f), hi(257, 257, 1.0f);
    FbmSettings s;
    s.octaves = 3; // niske oktave: rezolucija ih uhvati u potpunosti
    generateFbmTerrain(lo, s, 1.0f);
    generateFbmTerrain(hi, s, 1.0f);
    double maxDiff = 0.0;
    for (int y = 0; y < 129; y += 8)
        for (int x = 0; x < 129; x += 8)
            maxDiff = std::max(maxDiff, static_cast<double>(std::abs(lo.b[lo.idx(x, y)] - hi.b[hi.idx(2 * x, 2 * y)])));
    // Nije egzaktno: normalizacija koristi min/max mreže, a gušća mreža
    // uhvati malo ekstremnije vrijednosti. Tražimo slaganje na 2% opsega.
    CHECK(maxDiff < 0.02);
}

TEST(gradient_noise_je_nula_na_cjelobrojnim_tackama)
{
    // Osobina gradient šuma: u čvorovima rešetke vrijednost je 0.
    for (int i = -3; i <= 3; ++i)
        CHECK_NEAR(gradientNoise(static_cast<float>(i), static_cast<float>(2 * i), 42), 0.0f, 1e-6);
}

TEST(pgm_16bit_povratno_putovanje)
{
    Grid g(64, 48, 1.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);

    const std::string path = tempPath("prs_test_16.pgm");
    CHECK(writePGM(path, g.b.data(), g.width, g.height, 0.0f, 80.0f, true));

    const Heightmap hm = readPGM(path);
    CHECK(hm.width == 64);
    CHECK(hm.height == 48);
    float maxErr = 0.0f;
    for (std::size_t i = 0; i < g.cellCount(); ++i)
        maxErr = std::max(maxErr, std::abs(hm.values[i] * 80.0f - g.b[i]));
    CHECK(maxErr <= 80.0f / 65535.0f); // greška najviše jedan kvantizacijski korak
    std::filesystem::remove(path);
}

TEST(pgm_8bit_povratno_putovanje)
{
    const float data[6] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f, 2.0f}; // 2.0 se odsijeca na 1.0
    const std::string path = tempPath("prs_test_8.pgm");
    CHECK(writePGM(path, data, 3, 2, 0.0f, 1.0f, false));
    const Heightmap hm = readPGM(path);
    CHECK(hm.width == 3 && hm.height == 2);
    for (int i = 0; i < 5; ++i)
        CHECK_NEAR(hm.values[i], data[i], 0.5 / 255.0 + 1e-6);
    CHECK_NEAR(hm.values[5], 1.0f, 1e-6);
    std::filesystem::remove(path);
}

TEST(pgm_cita_zaglavlje_s_komentarom)
{
    const std::string path = tempPath("prs_test_comment.pgm");
    {
        std::FILE* f = std::fopen(path.c_str(), "wb");
        CHECK(f != nullptr);
        if (!f) return;
        const char header[] = "P5\n# komentar\n2 1\n255\n";
        std::fwrite(header, 1, sizeof(header) - 1, f);
        const unsigned char px[2] = {0, 255};
        std::fwrite(px, 1, 2, f);
        std::fclose(f);
    }
    const Heightmap hm = readPGM(path);
    CHECK(hm.width == 2 && hm.height == 1);
    CHECK_NEAR(hm.values[0], 0.0f, 1e-6);
    CHECK_NEAR(hm.values[1], 1.0f, 1e-6);
    std::filesystem::remove(path);
}

TEST(set_terrain_from_normalized_skalira)
{
    Grid g(2, 2, 1.0f);
    const float n[4] = {0.0f, 0.5f, 1.0f, 0.25f};
    setTerrainFromNormalized(g, n, 80.0f);
    CHECK_NEAR(g.b[1], 40.0f, 1e-6);
    CHECK_NEAR(g.b[2], 80.0f, 1e-6);
}

TEST(dijagnostika_mase)
{
    Grid g(10, 10, 2.0f);
    fillConstant(g.b, 1.0f);
    fillConstant(g.d, 0.5f);
    fillConstant(g.s, 0.25f);
    CHECK_NEAR(totalWater(g), 100 * 0.5 * 4.0, 1e-9);           // 100 ćelija × 0.5 × površina 4
    CHECK_NEAR(totalMaterial(g), 100 * (1.0 + 0.25) * 4.0, 1e-9);
}
