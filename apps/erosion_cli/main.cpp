// erosion_cli: generiše ili učitava teren i upisuje ga kao PGM.
// U koraku 2+ ovdje dodajemo pokretanje simulacije (--steps).

#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/ImageIO.h"
#include "erosion/Params.h"
#include "erosion/TerrainGen.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>
#include <string_view>

namespace {

void printUsage()
{
    std::puts(
        "Upotreba: erosion_cli [opcije]\n"
        "  --size N          mreža N x N (zadano 256)\n"
        "  --width W         širina mreže (nadjačava --size)\n"
        "  --height H        visina mreže (nadjačava --size)\n"
        "  --seed S          seed generatora (zadano 42)\n"
        "  --freq F          osnovna frekvencija šuma (zadano 4)\n"
        "  --octaves O       broj oktava fBm (zadano 6)\n"
        "  --persistence P   (zadano 0.5)\n"
        "  --lacunarity L    (zadano 2.0)\n"
        "  --exponent E      h = h^E nakon normalizacije (zadano 1.0)\n"
        "  --height-scale H  maksimalna visina terena (zadano 80)\n"
        "  --cell-size C     razmak ćelija (zadano 2)\n"
        "  --input file.pgm  učitaj teren iz PGM-a umjesto generisanja\n"
        "  --out file.pgm    izlazna heightmapa (zadano terrain.pgm)\n"
        "  --8bit            8-bitni PGM (zadano 16-bitni)\n"
        "  --help");
}

} // namespace

int main(int argc, char** argv)
{
    using namespace erosion;

    int         width  = 256;
    int         height = 256;
    FbmSettings fbmSettings;
    Params      params;
    std::string inputPath;
    std::string outPath = "terrain.pgm";
    bool        bits16  = true;

    try {
        for (int i = 1; i < argc; ++i) {
            const std::string_view arg = argv[i];
            auto next = [&]() -> const char* {
                if (i + 1 >= argc)
                    throw std::invalid_argument("nedostaje vrijednost za " + std::string(arg));
                return argv[++i];
            };

            if (arg == "--help" || arg == "-h") { printUsage(); return 0; }
            else if (arg == "--size")         { width = height = std::stoi(next()); }
            else if (arg == "--width")        { width  = std::stoi(next()); }
            else if (arg == "--height")       { height = std::stoi(next()); }
            else if (arg == "--seed")         { fbmSettings.seed = static_cast<uint32_t>(std::stoul(next())); }
            else if (arg == "--freq")         { fbmSettings.baseFreq    = std::stof(next()); }
            else if (arg == "--octaves")      { fbmSettings.octaves     = std::stoi(next()); }
            else if (arg == "--persistence")  { fbmSettings.persistence = std::stof(next()); }
            else if (arg == "--lacunarity")   { fbmSettings.lacunarity  = std::stof(next()); }
            else if (arg == "--exponent")     { fbmSettings.exponent    = std::stof(next()); }
            else if (arg == "--height-scale") { params.heightScale = std::stof(next()); }
            else if (arg == "--cell-size")    { params.cellSize    = std::stof(next()); }
            else if (arg == "--input")        { inputPath = next(); }
            else if (arg == "--out")          { outPath   = next(); }
            else if (arg == "--8bit")         { bits16 = false; }
            else {
                std::fprintf(stderr, "Nepoznata opcija: %s\n\n", argv[i]);
                printUsage();
                return 2;
            }
        }

        const auto t0 = std::chrono::steady_clock::now();

        Grid grid;
        if (!inputPath.empty()) {
            const Heightmap hm = readPGM(inputPath);
            grid.resize(hm.width, hm.height, params.cellSize);
            setTerrainFromNormalized(grid, hm.values.data(), params.heightScale);
        } else {
            grid.resize(width, height, params.cellSize);
            generateFbmTerrain(grid, fbmSettings, params.heightScale);
        }

        const auto   t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        const FieldRange r = fieldRange(grid.b);
        std::printf("Mreža:        %d x %d (%zu ćelija), cellSize = %.3f\n",
                    grid.width, grid.height, grid.cellCount(), grid.cellSize);
        std::printf("Teren:        %s\n", inputPath.empty() ? "fBm" : inputPath.c_str());
        std::printf("Visine:       min = %.3f, max = %.3f\n", r.min, r.max);
        std::printf("Materijal:    %.6e\n", totalMaterial(grid));
        std::printf("Vrijeme:      %.2f ms\n", ms);

        // Fiksni opseg [0, heightScale], da se slike prije i poslije erozije mogu porediti.
        if (!writePGM(outPath, grid.b.data(), grid.width, grid.height, 0.0f, params.heightScale, bits16)) {
            std::fprintf(stderr, "Greška: ne mogu upisati %s\n", outPath.c_str());
            return 1;
        }
        std::printf("Upisano:      %s\n", outPath.c_str());
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Greška: %s\n", e.what());
        return 1;
    }
    return 0;
}
