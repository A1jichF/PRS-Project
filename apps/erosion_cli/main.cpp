// erosion_cli: pravi početno stanje (teren + voda), pokreće simulaciju,
// ispisuje mjerenja i upisuje rezultate kao PGM.

#include "erosion/Bedrock.h"
#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/ImageIO.h"
#include "erosion/Params.h"
#include "erosion/Scenarios.h"
#include "erosion/Simulation.h"
#include "erosion/TerrainGen.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

void printUsage()
{
    std::puts(
        "Upotreba: erosion_cli [opcije]\n"
        "\n"
        "Mreža i teren:\n"
        "  --size N            mreža N x N (zadano 256)\n"
        "  --width W           širina mreže (nadjačava --size)\n"
        "  --height H          visina mreže (nadjačava --size)\n"
        "  --cell-size C       razmak ćelija (zadano 2)\n"
        "  --height-scale H    maksimalna visina terena (zadano 80)\n"
        "  --seed S            seed generatora (zadano 42)\n"
        "  --freq F            osnovna frekvencija šuma (zadano 4)\n"
        "  --octaves O         broj oktava fBm (zadano 6)\n"
        "  --persistence P     (zadano 0.5)\n"
        "  --lacunarity L      (zadano 2.0)\n"
        "  --exponent E        h = h^E nakon normalizacije (zadano 1.0)\n"
        "  --input file.pgm    učitaj teren iz PGM-a umjesto generisanja\n"
        "\n"
        "Početno stanje vode (--scenario):\n"
        "  fbm                 bez vode (zadano), voda dolazi od kiše/izvora\n"
        "  layer               sloj vode --depth svuda preko terena\n"
        "  lake                voda do nivoa --level (mirno jezero)\n"
        "  dambreak            ravan teren, voda --depth i stub 3*depth u centru\n"
        "  river               dolina s vijugavim koritom od stijene, izvorom na vrhu\n"
        "                      i bazenom na dnu (teren se ne generiše iz fBm-a)\n"
        "  --depth D           dubina za layer/dambreak (zadano 1)\n"
        "  --level L           nivo vode za lake (zadano heightScale/2)\n"
        "\n"
        "Simulacija:\n"
        "  --steps N           broj koraka (zadano 0: samo početno stanje)\n"
        "  --dt T              vremenski korak (zadano 0.02)\n"
        "  --rain R            kiša: visina vode po ćeliji u sekundi (zadano 0)\n"
        "  --evaporation K     udio vode koji ispari u sekundi (zadano 0.015, 0 = bez)\n"
        "  --source X Y R Q    izvor vode u (X, Y), poluprečnik R, Q po sekundi\n"
        "                      (može se ponoviti)\n"
        "\n"
        "Stijena (bedrock); zadaje se poslije terena, redom kako je navedeno:\n"
        "  --soil-depth D      isti sloj tla debljine D svuda iznad stijene\n"
        "                      (za river: debljina tla izvan korita, zadano 1.5)\n"
        "  --rock-level L      stijena do visine L (ispod L je gola stijena)\n"
        "  --bedrock file.pgm  mapa visina stijene (skalirana s --height-scale)\n"
        "  --no-bedrock        ukloni stijenu (npr. river bez stijene, za poređenje)\n"
        "  --rock-erodibility K  tvrdoća: 0 = stijena se ne erodira (zadano),\n"
        "                      0.05 = 20x sporije od tla, 1 = kao tlo\n"
        "  --river-rate Q      jačina izvora za river (zadano 0.6)\n"
        "\n"
        "Erozija:\n"
        "  --no-erosion        isključi eroziju i depoziciju\n"
        "  --kc K              kapacitet nošenja sedimenta (zadano 0.5)\n"
        "  --ks K              brzina erozije [1/s] (zadano 0.5)\n"
        "  --kd K              brzina depozicije [1/s] (zadano 1.0)\n"
        "  --kdmax D           dubina do koje raste kapacitet, Šťava 2008 (zadano 0.5)\n"
        "  --transport T       transport sedimenta: flux (zadano, čuva masu)\n"
        "                      ili sl (Mei et al., semi-Lagrangian, gubi masu)\n"
        "\n"
        "Izlaz:\n"
        "  --out file.pgm      teren (zadano terrain.pgm)\n"
        "  --water-out f.pgm   dubina vode (zadano water.pgm)\n"
        "  --diff-out f.pgm    promjena terena b - b0 (sivo = bez promjene,\n"
        "                      tamno = erozija, svijetlo = taloženje)\n"
        "  --sediment-out f    sediment u vodi\n"
        "  --soil-out f.pgm    debljina tla iznad stijene (crno = gola stijena)\n"
        "  --snapshot-every K  upiši vodu svakih K koraka kao <prefix>_NNNNNN.pgm\n"
        "  --snapshot-prefix P (zadano water)\n"
        "  --8bit              8-bitni PGM (zadano 16-bitni)\n"
        "  --csv file.csv      dopiši red s mjerenjem (ms/korak po fazama); zaglavlje\n"
        "                      se upisuje ako fajl ne postoji\n"
        "  --tag T             oznaka u CSV redu (npr. seq, omp8), zadano seq\n"
        "  --help");
}

enum class Scenario { Fbm, Layer, Lake, DamBreak, River };

Scenario parseScenario(std::string_view s)
{
    if (s == "fbm")      return Scenario::Fbm;
    if (s == "layer")    return Scenario::Layer;
    if (s == "lake")     return Scenario::Lake;
    if (s == "dambreak") return Scenario::DamBreak;
    if (s == "river")    return Scenario::River;
    throw std::invalid_argument("nepoznat scenario: " + std::string(s));
}

} // namespace

int main(int argc, char** argv)
{
    using namespace erosion;

    int           width  = 256;
    int           height = 256;
    FbmSettings   fbmSettings;
    WorldSettings world;
    Params        params;
    std::vector<WaterSource> sources;

    Scenario    scenario = Scenario::Fbm;
    float       depth    = 1.0f;
    float       level    = -1.0f; // < 0 znači heightScale / 2
    int         steps    = 0;
    std::string inputPath;
    std::string outPath        = "terrain.pgm";
    std::string waterOutPath   = "water.pgm";
    std::string diffOutPath;
    std::string sedimentOutPath;
    std::string soilOutPath;
    std::string bedrockPath;
    RiverSettings river;
    const float kUnset     = std::numeric_limits<float>::quiet_NaN();
    float       soilDepth  = kUnset;
    float       rockLevel  = kUnset;
    bool        noBedrock  = false;
    std::string csvPath;
    std::string csvTag = "seq";
    std::string snapshotPrefix = "water";
    int         snapshotEvery  = 0;
    bool        bits16         = true;

    try {
        for (int i = 1; i < argc; ++i) {
            const std::string_view arg = argv[i];
            auto next = [&]() -> const char* {
                if (i + 1 >= argc)
                    throw std::invalid_argument("nedostaje vrijednost za " + std::string(arg));
                return argv[++i];
            };

            if (arg == "--help" || arg == "-h") { printUsage(); return 0; }
            else if (arg == "--size")            { width = height = std::stoi(next()); }
            else if (arg == "--width")           { width  = std::stoi(next()); }
            else if (arg == "--height")          { height = std::stoi(next()); }
            else if (arg == "--cell-size")       { world.cellSize    = std::stof(next()); }
            else if (arg == "--height-scale")    { world.heightScale = std::stof(next()); }
            else if (arg == "--seed")            { fbmSettings.seed = static_cast<uint32_t>(std::stoul(next())); }
            else if (arg == "--freq")            { fbmSettings.baseFreq    = std::stof(next()); }
            else if (arg == "--octaves")         { fbmSettings.octaves     = std::stoi(next()); }
            else if (arg == "--persistence")     { fbmSettings.persistence = std::stof(next()); }
            else if (arg == "--lacunarity")      { fbmSettings.lacunarity  = std::stof(next()); }
            else if (arg == "--exponent")        { fbmSettings.exponent    = std::stof(next()); }
            else if (arg == "--input")           { inputPath = next(); }
            else if (arg == "--scenario")        { scenario = parseScenario(next()); }
            else if (arg == "--depth")           { depth = std::stof(next()); }
            else if (arg == "--level")           { level = std::stof(next()); }
            else if (arg == "--steps")           { steps = std::stoi(next()); }
            else if (arg == "--dt")              { params.dt = std::stof(next()); }
            else if (arg == "--rain")            { params.rainRate = std::stof(next()); }
            else if (arg == "--evaporation")     { params.evaporationRate = std::stof(next()); }
            else if (arg == "--source") {
                WaterSource s;
                s.x      = std::stoi(next());
                s.y      = std::stoi(next());
                s.radius = std::stof(next());
                s.rate   = std::stof(next());
                sources.push_back(s);
            }
            else if (arg == "--soil-depth")      { soilDepth = std::stof(next()); }
            else if (arg == "--rock-level")      { rockLevel = std::stof(next()); }
            else if (arg == "--bedrock")         { bedrockPath = next(); }
            else if (arg == "--no-bedrock")      { noBedrock = true; }
            else if (arg == "--rock-erodibility") { params.rockErodibility = std::stof(next()); }
            else if (arg == "--river-rate")      { river.sourceRate = std::stof(next()); }
            else if (arg == "--no-erosion")      { params.erosionEnabled   = false; }
            else if (arg == "--kc")              { params.sedimentCapacity = std::stof(next()); }
            else if (arg == "--ks")              { params.dissolvingRate   = std::stof(next()); }
            else if (arg == "--kd")              { params.depositionRate   = std::stof(next()); }
            else if (arg == "--kdmax")           { params.maxErosionDepth  = std::stof(next()); }
            else if (arg == "--transport") {
                const std::string_view t = next();
                if (t == "flux")    params.sedimentTransport = SedimentTransport::Flux;
                else if (t == "sl") params.sedimentTransport = SedimentTransport::SemiLagrangian;
                else throw std::invalid_argument("nepoznat transport: " + std::string(t) + " (flux ili sl)");
            }
            else if (arg == "--out")             { outPath = next(); }
            else if (arg == "--water-out")       { waterOutPath = next(); }
            else if (arg == "--diff-out")        { diffOutPath = next(); }
            else if (arg == "--sediment-out")    { sedimentOutPath = next(); }
            else if (arg == "--soil-out")        { soilOutPath = next(); }
            else if (arg == "--snapshot-every")  { snapshotEvery = std::stoi(next()); }
            else if (arg == "--snapshot-prefix") { snapshotPrefix = next(); }
            else if (arg == "--8bit")            { bits16 = false; }
            else if (arg == "--csv")             { csvPath = next(); }
            else if (arg == "--tag")             { csvTag = next(); }
            else {
                std::fprintf(stderr, "Nepoznata opcija: %s\n\n", argv[i]);
                printUsage();
                return 2;
            }
        }

        // ── Početno stanje ───────────────────────────────────────────────────
        Grid grid;
        if (!inputPath.empty()) {
            const Heightmap hm = readPGM(inputPath);
            grid.resize(hm.width, hm.height, world.cellSize);
            setTerrainFromNormalized(grid, hm.values.data(), world.heightScale);
        } else {
            grid.resize(width, height, world.cellSize);
            if (scenario != Scenario::DamBreak && scenario != Scenario::River)
                generateFbmTerrain(grid, fbmSettings, world.heightScale);
        }

        switch (scenario) {
        case Scenario::Fbm:      break;
        case Scenario::Layer:    fillConstant(grid.d, depth); break;
        case Scenario::Lake:     fillToLevel(grid, level >= 0.0f ? level : 0.5f * world.heightScale); break;
        case Scenario::DamBreak: setupDamBreak(grid, depth, 3.0f * depth, std::max(1, std::min(grid.width, grid.height) / 8)); break;
        case Scenario::River:
            // Teren, stijena (korito) i izvor dolaze iz scenarija.
            if (!std::isnan(soilDepth))
                river.soilDepth = soilDepth;
            sources.push_back(setupRiverValley(grid, river));
            break;
        }

        // ── Stijena ──────────────────────────────────────────────────────────
        if (scenario != Scenario::River && !std::isnan(soilDepth))
            setBedrockUniformSoil(grid, soilDepth);
        if (!std::isnan(rockLevel))
            setBedrockLevel(grid, rockLevel);
        if (!bedrockPath.empty()) {
            const Heightmap hm = readPGM(bedrockPath);
            if (hm.width != grid.width || hm.height != grid.height)
                throw std::invalid_argument("--bedrock: dimenzije mape ne odgovaraju mreži");
            std::vector<float> heights(hm.values.size());
            for (std::size_t i = 0; i < heights.size(); ++i)
                heights[i] = hm.values[i] * world.heightScale;
            setBedrockHeights(grid, heights.data());
        }
        if (noBedrock)
            clearBedrock(grid);
        const BedrockStats rock0 = bedrockStats(grid);

        const double water0    = totalWater(grid);
        const double material0 = totalMaterial(grid);
        const Field  b0        = grid.b;

        Simulation sim(std::move(grid), params);
        for (const WaterSource& s : sources)
            sim.addSource(s);

        std::printf("Mreža:          %d x %d (%zu ćelija), cellSize = %.3f\n",
                    sim.grid().width, sim.grid().height, sim.grid().cellCount(), sim.grid().cellSize);
        std::printf("Koraci:         %d x dt = %.4f s\n", steps, static_cast<double>(params.dt));

        // ── Simulacija ───────────────────────────────────────────────────────
        float maxCfl = 0.0f;
        auto  writeSnapshot = [&](long stepIndex) {
            char name[512];
            std::snprintf(name, sizeof(name), "%s_%06ld.pgm", snapshotPrefix.c_str(), stepIndex);
            writePGM(name, sim.grid().d.data(), sim.grid().width, sim.grid().height, bits16);
        };
        if (snapshotEvery > 0)
            writeSnapshot(0);

        const auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < steps; ++i) {
            sim.step();
            if (snapshotEvery > 0 && (i + 1) % snapshotEvery == 0) {
                // CFL i snimci se računaju samo povremeno, da ne kvare mjerenje vremena.
                maxCfl = std::max(maxCfl, velocityCfl(sim.grid(), sim.params()).cfl);
                writeSnapshot(i + 1);
            }
        }
        const auto   t1    = std::chrono::steady_clock::now();
        const double wallMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
        maxCfl = std::max(maxCfl, velocityCfl(sim.grid(), sim.params()).cfl);

        // ── Izvještaj ────────────────────────────────────────────────────────
        const Grid&        g  = sim.grid();
        const PhaseTimings& t = sim.timings();
        if (steps > 0) {
            const double msPerStep = t.total() / static_cast<double>(t.steps);
            const double mcells    = static_cast<double>(g.cellCount()) / (msPerStep * 1e3);
            std::printf("Ukupno vrijeme: %.2f ms (%.4f ms/korak, %.1f Mćelija/s)\n", wallMs, msPerStep, mcells);
            std::printf("  dodavanje vode %.4f ms/korak\n", t.addWater / static_cast<double>(t.steps));
            std::printf("  fluks          %.4f ms/korak\n", t.flux / static_cast<double>(t.steps));
            std::printf("  voda + brzina  %.4f ms/korak\n", t.water / static_cast<double>(t.steps));
            std::printf("  erozija        %.4f ms/korak\n", t.erosion / static_cast<double>(t.steps));
            std::printf("  transport sed. %.4f ms/korak (%s)\n", t.advect / static_cast<double>(t.steps),
                        params.sedimentTransport == SedimentTransport::Flux ? "flux" : "semi-Lagrangian");
            std::printf("  isparavanje    %.4f ms/korak\n", t.evaporate / static_cast<double>(t.steps));

            if (!csvPath.empty()) {
                bool needHeader = true;
                if (std::FILE* probe = std::fopen(csvPath.c_str(), "rb")) {
                    needHeader = std::fgetc(probe) == EOF;
                    std::fclose(probe);
                }
                std::FILE* csv = std::fopen(csvPath.c_str(), "ab");
                if (!csv)
                    throw std::runtime_error("ne mogu otvoriti " + csvPath);
                if (needHeader)
                    std::fputs("tag,width,height,steps,ms_per_step,mcells_per_s,"
                               "add_water,flux,transport,water,erosion,evaporate\n", csv);
                const double n = static_cast<double>(t.steps);
                std::fprintf(csv, "%s,%d,%d,%ld,%.6f,%.3f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                             csvTag.c_str(), g.width, g.height, t.steps, msPerStep, mcells,
                             t.addWater / n, t.flux / n, t.advect / n, t.water / n,
                             t.erosion / n, t.evaporate / n);
                std::fclose(csv);
            }
        }

        const FieldRange dr = fieldRange(g.d);
        const double     water1 = totalWater(g);
        std::printf("Voda:           %.6e -> %.6e", water0, water1);
        if (water0 > 0.0 && params.rainRate == 0.0f && sources.empty()) {
            if (params.evaporationRate == 0.0f) {
                std::printf("  (relativna promjena %.2e)", (water1 - water0) / water0);
            } else {
                // Bez kiše i izvora voda se smanjuje tačno faktorom isparavanja.
                const double expected =
                    water0 * std::pow(std::max(0.0, 1.0 - double(params.evaporationRate) * params.dt), steps);
                std::printf("  (očekivano zbog isparavanja %.6e, odstupanje %.2e)",
                            expected, expected > 0.0 ? (water1 - expected) / expected : 0.0);
            }
        }
        std::printf("\n");
        std::printf("Dubina vode:    min = %.6f, max = %.6f\n", dr.min, dr.max);
        const double material1 = totalMaterial(g);
        std::printf("Materijal:      %.6e -> %.6e", material0, material1);
        if (material0 > 0.0)
            std::printf("  (relativna promjena %.2e)", (material1 - material0) / material0);
        std::printf("\n");

        if (rock0.rockCells > 0) {
            const BedrockStats rock1 = bedrockStats(g);
            const double       cells = static_cast<double>(g.cellCount());
            std::printf("Stijena:        zadana u %.1f %% ćelija, tvrdoća %.3f\n",
                        100.0 * static_cast<double>(rock0.rockCells) / cells,
                        static_cast<double>(params.rockErodibility));
            std::printf("  gola stijena  %.2f %% -> %.2f %% ćelija\n",
                        100.0 * static_cast<double>(rock0.bareCells) / cells,
                        100.0 * static_cast<double>(rock1.bareCells) / cells);
            std::printf("  tlo           %.4e -> %.4e (zapremina)\n", rock0.soilVolume, rock1.soilVolume);
        }

        // Erozija: koliko je terena odneseno, a koliko nataloženo.
        Field  diff(g.cellCount());
        double eroded = 0.0, deposited = 0.0;
        for (std::size_t i = 0; i < g.cellCount(); ++i) {
            diff[i] = g.b[i] - b0[i];
            if (diff[i] < 0.0f) eroded    -= diff[i];
            else                deposited += diff[i];
        }
        const double area = static_cast<double>(g.cellSize) * g.cellSize;
        const FieldRange diffRange = fieldRange(diff);
        std::printf("Erozija:        odneseno %.4e, nataloženo %.4e, u vodi %.4e (zapremina)\n",
                    eroded * area, deposited * area, fieldSum(g.s) * area);
        std::printf("Promjena b:     min = %.4f, max = %.4f\n", diffRange.min, diffRange.max);
        std::printf("CFL (max):      %.4f%s\n", maxCfl,
                    maxCfl >= 0.999f ? "  (neke ćelije se isprazne u jednom koraku: tanak sloj na strmini)" : "");
        if (!allFinite(g.d) || !allFinite(g.u) || !allFinite(g.v))
            std::printf("UPOZORENJE: NaN ili beskonačne vrijednosti u stanju!\n");

        // ── Izlazne slike ────────────────────────────────────────────────────
        // Teren u fiksnom opsegu [0, heightScale], da se slike prije i poslije erozije mogu porediti.
        if (!writePGM(outPath, g.b.data(), g.width, g.height, 0.0f, world.heightScale, bits16))
            throw std::runtime_error("ne mogu upisati " + outPath);
        if (!writePGM(waterOutPath, g.d.data(), g.width, g.height, bits16))
            throw std::runtime_error("ne mogu upisati " + waterOutPath);
        std::printf("Upisano:        %s, %s", outPath.c_str(), waterOutPath.c_str());
        if (!diffOutPath.empty()) {
            // Simetričan opseg oko nule, da siva boja znači "bez promjene".
            const float m = std::max(std::abs(diffRange.min), std::abs(diffRange.max));
            const float r = m > 0.0f ? m : 1.0f;
            if (!writePGM(diffOutPath, diff.data(), g.width, g.height, -r, r, bits16))
                throw std::runtime_error("ne mogu upisati " + diffOutPath);
            std::printf(", %s", diffOutPath.c_str());
        }
        if (!soilOutPath.empty()) {
            // Debljina tla b - r; bez stijene je neograničena, pa se prikazuje kao najveća.
            Field soil(g.cellCount(), 0.0f);
            float maxSoil = 0.0f;
            for (std::size_t i = 0; i < g.cellCount(); ++i)
                if (g.r[i] != kNoBedrock) {
                    soil[i] = std::max(0.0f, g.b[i] - g.r[i]);
                    maxSoil = std::max(maxSoil, soil[i]);
                }
            for (std::size_t i = 0; i < g.cellCount(); ++i)
                if (g.r[i] == kNoBedrock)
                    soil[i] = maxSoil;
            if (!writePGM(soilOutPath, soil.data(), g.width, g.height, 0.0f,
                          maxSoil > 0.0f ? maxSoil : 1.0f, bits16))
                throw std::runtime_error("ne mogu upisati " + soilOutPath);
            std::printf(", %s", soilOutPath.c_str());
        }
        if (!sedimentOutPath.empty()) {
            if (!writePGM(sedimentOutPath, g.s.data(), g.width, g.height, bits16))
                throw std::runtime_error("ne mogu upisati " + sedimentOutPath);
            std::printf(", %s", sedimentOutPath.c_str());
        }
        std::printf("\n");
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Greška: %s\n", e.what());
        return 1;
    }
    return 0;
}
