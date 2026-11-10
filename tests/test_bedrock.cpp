// Testovi za stijenu (bedrock): polje Grid::r i sloj tla iznad njega.

#include "TestFramework.h"

#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/Kernels.h"
#include "erosion/Params.h"
#include "erosion/Simulation.h"
#include "erosion/TerrainGen.h"

#include <cstdint>
#include <iterator>

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
    // Dok kernel erozije ne koristi stijenu, r ostaje netaknuto kroz simulaciju.
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
