// C interfejs (prs_erosion.h) oko erosion::Simulation, za DLL.
//
// Pravila na granici DLL-a:
//   - nijedan C++ izuzetak ne smije izaći van (svaka funkcija koja alocira
//     ima try/catch i grešku vraća kroz prs_erosion_last_error);
//   - NULL handle je dozvoljen svuda i ne radi ništa (ili vraća 0 / NULL).

#include "prs_erosion.h"

#include "erosion/Diagnostics.h"
#include "erosion/Grid.h"
#include "erosion/Params.h"
#include "erosion/Simulation.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <memory>
#include <string>

struct PrsErosion {
    std::unique_ptr<erosion::Simulation> sim;
    double lastStepMs = 0.0;
};

namespace {

thread_local std::string tLastError;

void setError(std::string msg) { tLastError = std::move(msg); }
void clearError() { tLastError.clear(); }

PrsErosionParams toC(const erosion::Params& p)
{
    PrsErosionParams c{};
    c.dt                  = p.dt;
    c.gravity             = p.gravity;
    c.pipeArea            = p.pipeArea;
    c.rainRate            = p.rainRate;
    c.evaporationRate     = p.evaporationRate;
    c.minDepthForVelocity = p.minDepthForVelocity;
    c.erosionEnabled      = p.erosionEnabled ? 1 : 0;
    c.sedimentCapacity    = p.sedimentCapacity;
    c.dissolvingRate      = p.dissolvingRate;
    c.depositionRate      = p.depositionRate;
    c.minTilt             = p.minTilt;
    c.maxErosionDepth     = p.maxErosionDepth;
    c.sedimentTransport   = p.sedimentTransport == erosion::SedimentTransport::Flux
                                ? PRS_TRANSPORT_FLUX
                                : PRS_TRANSPORT_SEMI_LAGRANGIAN;
    return c;
}

/// Prevodi C parametre u erosion::Params. Vraća false (i postavlja grešku)
/// ako vrijednost nema smisla, da simulacija ne bi proizvela NaN.
bool fromC(const PrsErosionParams& c, erosion::Params& out)
{
    auto finite = [](float v) { return std::isfinite(v); };
    if (!finite(c.dt) || c.dt <= 0.0f) {
        setError("dt mora biti > 0");
        return false;
    }
    const float nonNegative[] = {c.gravity, c.pipeArea, c.rainRate, c.evaporationRate,
                                 c.minDepthForVelocity, c.sedimentCapacity,
                                 c.dissolvingRate, c.depositionRate, c.minTilt,
                                 c.maxErosionDepth};
    for (float v : nonNegative) {
        if (!finite(v) || v < 0.0f) {
            setError("parametri moraju biti konačni i >= 0");
            return false;
        }
    }
    if (c.sedimentTransport != PRS_TRANSPORT_FLUX &&
        c.sedimentTransport != PRS_TRANSPORT_SEMI_LAGRANGIAN) {
        setError("nepoznat sedimentTransport");
        return false;
    }

    erosion::Params p = out; // polja koja C API ne izlaže ostaju kakva jesu
    p.dt                  = c.dt;
    p.gravity             = c.gravity;
    p.pipeArea            = c.pipeArea;
    p.rainRate            = c.rainRate;
    p.evaporationRate     = c.evaporationRate;
    p.minDepthForVelocity = c.minDepthForVelocity;
    p.erosionEnabled      = c.erosionEnabled != 0;
    p.sedimentCapacity    = c.sedimentCapacity;
    p.dissolvingRate      = c.dissolvingRate;
    p.depositionRate      = c.depositionRate;
    p.minTilt             = c.minTilt;
    p.maxErosionDepth     = c.maxErosionDepth;
    p.sedimentTransport   = c.sedimentTransport == PRS_TRANSPORT_FLUX
                                ? erosion::SedimentTransport::Flux
                                : erosion::SedimentTransport::SemiLagrangian;
    out = p;
    return true;
}

erosion::Grid makeGrid(int width, int height, float cellSize, const float* terrain)
{
    erosion::Grid g(width, height, cellSize);
    if (terrain)
        std::copy(terrain, terrain + g.cellCount(), g.b.begin());
    return g;
}

} // namespace

extern "C" {

int prs_erosion_api_version(void) { return PRS_EROSION_API_VERSION; }

const char* prs_erosion_last_error(void) { return tLastError.c_str(); }

void prs_erosion_default_params(PrsErosionParams* out)
{
    if (out)
        *out = toC(erosion::Params{});
}

PrsErosion* prs_erosion_create(int width, int height, float cellSize, const float* terrain,
                               const PrsErosionParams* params)
{
    clearError();
    if (width < 2 || height < 2) {
        setError("mreža mora biti barem 2 x 2");
        return nullptr;
    }
    if (!std::isfinite(cellSize) || cellSize <= 0.0f) {
        setError("cellSize mora biti > 0");
        return nullptr;
    }
    erosion::Params p;
    if (params && !fromC(*params, p))
        return nullptr;

    try {
        auto h = std::make_unique<PrsErosion>();
        h->sim = std::make_unique<erosion::Simulation>(
            makeGrid(width, height, cellSize, terrain), p);
        return h.release();
    } catch (const std::exception& e) {
        setError(std::string("prs_erosion_create: ") + e.what());
    } catch (...) {
        setError("prs_erosion_create: nepoznata greška");
    }
    return nullptr;
}

void prs_erosion_destroy(PrsErosion* sim) { delete sim; }

int prs_erosion_reset(PrsErosion* sim, const float* terrain)
{
    clearError();
    if (!sim) {
        setError("NULL handle");
        return 0;
    }
    try {
        const erosion::Grid& old = sim->sim->grid();
        auto fresh = std::make_unique<erosion::Simulation>(
            makeGrid(old.width, old.height, old.cellSize, terrain), sim->sim->params());
        for (const erosion::WaterSource& src : sim->sim->sources())
            fresh->addSource(src);
        sim->sim        = std::move(fresh);
        sim->lastStepMs = 0.0;
        return 1;
    } catch (const std::exception& e) {
        setError(std::string("prs_erosion_reset: ") + e.what());
    } catch (...) {
        setError("prs_erosion_reset: nepoznata greška");
    }
    return 0;
}

int prs_erosion_set_params(PrsErosion* sim, const PrsErosionParams* params)
{
    clearError();
    if (!sim || !params) {
        setError("NULL argument");
        return 0;
    }
    return fromC(*params, sim->sim->params()) ? 1 : 0;
}

void prs_erosion_get_params(const PrsErosion* sim, PrsErosionParams* out)
{
    if (sim && out)
        *out = toC(sim->sim->params());
}

int prs_erosion_add_source(PrsErosion* sim, int x, int y, float radius, float rate)
{
    clearError();
    if (!sim) {
        setError("NULL handle");
        return 0;
    }
    if (!std::isfinite(radius) || radius < 0.0f || !std::isfinite(rate)) {
        setError("neispravan izvor vode");
        return 0;
    }
    try {
        sim->sim->addSource(erosion::WaterSource{x, y, radius, rate});
        return 1;
    } catch (...) {
        setError("prs_erosion_add_source: nema memorije");
    }
    return 0;
}

void prs_erosion_clear_sources(PrsErosion* sim)
{
    if (sim)
        sim->sim->clearSources();
}

void prs_erosion_step(PrsErosion* sim, int steps)
{
    if (!sim || steps <= 0)
        return;
    const auto t0 = std::chrono::steady_clock::now();
    sim->sim->run(steps);
    const auto t1 = std::chrono::steady_clock::now();
    sim->lastStepMs = std::chrono::duration<double, std::milli>(t1 - t0).count() / steps;
}

int prs_erosion_width(const PrsErosion* sim) { return sim ? sim->sim->grid().width : 0; }
int prs_erosion_height(const PrsErosion* sim) { return sim ? sim->sim->grid().height : 0; }
float prs_erosion_cell_size(const PrsErosion* sim) { return sim ? sim->sim->grid().cellSize : 0.0f; }

const float* prs_erosion_terrain(const PrsErosion* sim) { return sim ? sim->sim->grid().b.data() : nullptr; }
const float* prs_erosion_water(const PrsErosion* sim) { return sim ? sim->sim->grid().d.data() : nullptr; }
const float* prs_erosion_sediment(const PrsErosion* sim) { return sim ? sim->sim->grid().s.data() : nullptr; }
const float* prs_erosion_velocity_x(const PrsErosion* sim) { return sim ? sim->sim->grid().u.data() : nullptr; }
const float* prs_erosion_velocity_y(const PrsErosion* sim) { return sim ? sim->sim->grid().v.data() : nullptr; }

void prs_erosion_get_stats(const PrsErosion* sim, PrsErosionStats* out)
{
    if (!out)
        return;
    *out = PrsErosionStats{};
    if (!sim)
        return;
    const erosion::Grid& g = sim->sim->grid();
    out->steps         = sim->sim->stepCount();
    out->simTime       = sim->sim->time();
    out->lastStepMs    = sim->lastStepMs;
    out->totalWater    = erosion::totalWater(g);
    out->totalMaterial = erosion::totalMaterial(g);
}

} // extern "C"
