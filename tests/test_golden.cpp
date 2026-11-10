// Zaštita od regresije: kontrolni zbirovi polja nakon referentne simulacije.
//
// Vrijednosti su snimljene na commitu 0c710d6 (prije uvođenja stijene) i
// služe kao dokaz da nova mogućnost, dok se ne koristi, ne mijenja NI JEDAN
// BIT rezultata. Svaka izmjena kernela koja nenamjerno promijeni rezultat
// obori ovaj test.
//
// Zbirovi zavise od redoslijeda float operacija, pa važe za sekvencijalnu
// verziju prevedenu bez -ffast-math i bez FMA (-march=native). Uz te uslove
// su isti u Release (-O2) i Debug (-O0) buildu (provjereno, GCC 16.2). Ako se
// promijene opcije kompajlera ili se rezultat mijenja NAMJERNO, nove
// vrijednosti (test ih ispiše kad padne) treba upisati ovdje i to navesti u
// poruci commita.

#include "TestFramework.h"

#include "erosion/Grid.h"
#include "erosion/Kernels.h"
#include "erosion/Params.h"
#include "erosion/Simulation.h"
#include "erosion/TerrainGen.h"

#include <cstdint>
#include <cstdio>
#include <cstring>

using namespace erosion;

namespace {

/// FNV-1a (64 bita) preko sirovih bajtova polja: mijenja se ako se promijeni
/// ijedan bit ijedne vrijednosti.
std::uint64_t fieldHash(const Field& f)
{
    std::uint64_t h = 14695981039346656037ull;
    for (float value : f) {
        unsigned char bytes[sizeof(float)];
        std::memcpy(bytes, &value, sizeof bytes);
        for (unsigned char byte : bytes) {
            h ^= byte;
            h *= 1099511628211ull;
        }
    }
    return h;
}

struct GoldenHashes {
    std::uint64_t b, d, s;
};

/// fBm 96 x 96, kiša + jedan izvor, 300 koraka.
Simulation runReference(SedimentTransport transport)
{
    Grid g(96, 96, 2.0f);
    generateFbmTerrain(g, FbmSettings{}, 80.0f);

    Params p;
    p.rainRate          = 0.01f;
    p.sedimentTransport = transport;

    Simulation sim(std::move(g), p);
    sim.addSource(WaterSource{30, 40, 4.0f, 2.0f});
    sim.run(300);
    return sim;
}

void checkGolden(const char* name, const Simulation& sim, const GoldenHashes& expected)
{
    const Grid&        g = sim.grid();
    const GoldenHashes actual{fieldHash(g.b), fieldHash(g.d), fieldHash(g.s)};
    const bool ok = actual.b == expected.b && actual.d == expected.d && actual.s == expected.s;
    CHECK(ok);
    if (!ok)
        std::fprintf(stderr,
                     "  %s: dobijeno {0x%016llxull, 0x%016llxull, 0x%016llxull}\n", name,
                     static_cast<unsigned long long>(actual.b),
                     static_cast<unsigned long long>(actual.d),
                     static_cast<unsigned long long>(actual.s));
}

} // namespace

TEST(golden_flux_transport_bit_identican)
{
    checkGolden("flux", runReference(SedimentTransport::Flux),
                {0xe204dc3c55ed72beull, 0x6a7978d561a7617dull, 0x17ef2d188a375a22ull});
}

TEST(golden_semi_lagrangian_bit_identican)
{
    checkGolden("semi-lagrangian", runReference(SedimentTransport::SemiLagrangian),
                {0x54e047b145ef9529ull, 0x1e56aedce2d87feeull, 0x9770c6c8cfe3ecd0ull});
}
