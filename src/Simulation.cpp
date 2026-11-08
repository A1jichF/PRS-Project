#include "erosion/Simulation.h"

#include <chrono>
#include <utility>

namespace erosion {

namespace {

using Clock = std::chrono::steady_clock;

double msSince(Clock::time_point t0, Clock::time_point t1)
{
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

} // namespace

Simulation::Simulation(Grid grid, Params params)
    : m_grid(std::move(grid))
    , m_params(params)
{
    allocateScratch();
}

void Simulation::allocateScratch()
{
    // Pomoćni baferi u slotovima 9–11, različitim od polja mreže (0–8),
    // da se ne poklope s nizovima uz koje se koriste (vidi AlignedAllocator).
    StaggerScope stagger(9);
    for (Field* f : {&m_bScratch, &m_sScratch, &m_rateScratch}) {
        Field().swap(*f);
        f->assign(m_grid.cellCount(), 0.0f);
    }
}

void Simulation::step()
{
    // Ako je neko promijenio dimenzije mreže preko grid(), uskladi bafere.
    if (m_bScratch.size() != m_grid.cellCount())
        allocateScratch();

    const bool erosion = m_params.erosionEnabled;
    const bool fluxTransport =
        erosion && m_params.sedimentTransport == SedimentTransport::Flux;
    const bool slTransport =
        erosion && m_params.sedimentTransport == SedimentTransport::SemiLagrangian;

    const auto t0 = Clock::now();
    seq::addWater(m_grid, m_params, m_sources);
    const auto t1 = Clock::now();
    seq::computeFlux(m_grid, m_params);
    const auto t2 = Clock::now();
    // Konzervativni transport treba d1 (vodu prije toka), pa ide prije ažuriranja vode.
    if (fluxTransport)
        seq::transportSedimentFlux(m_grid, m_params, m_sScratch, m_rateScratch);
    const auto t3 = Clock::now();
    seq::updateWaterAndVelocity(m_grid, m_params);
    const auto t4 = Clock::now();
    if (erosion)
        seq::erodeDeposit(m_grid, m_params, m_bScratch);
    const auto t5 = Clock::now();
    // Mei et al.: advekcija ide na kraju, s brzinom iz ovog koraka.
    if (slTransport)
        seq::advectSediment(m_grid, m_params, m_sScratch);
    const auto t6 = Clock::now();
    seq::evaporate(m_grid, m_params);
    const auto t7 = Clock::now();

    m_timings.evaporate += msSince(t6, t7);
    m_timings.addWater += msSince(t0, t1);
    m_timings.flux     += msSince(t1, t2);
    m_timings.advect   += msSince(t2, t3) + msSince(t5, t6);
    m_timings.water    += msSince(t3, t4);
    m_timings.erosion  += msSince(t4, t5);
    ++m_timings.steps;

    ++m_stepCount;
    m_time += m_params.dt;
}

void Simulation::run(int n)
{
    for (int i = 0; i < n; ++i)
        step();
}

} // namespace erosion
