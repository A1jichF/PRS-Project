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
}

void Simulation::step()
{
    const auto t0 = Clock::now();
    seq::addWater(m_grid, m_params, m_sources);
    const auto t1 = Clock::now();
    seq::computeFlux(m_grid, m_params);
    const auto t2 = Clock::now();
    seq::updateWaterAndVelocity(m_grid, m_params);
    const auto t3 = Clock::now();

    m_timings.addWater += msSince(t0, t1);
    m_timings.flux     += msSince(t1, t2);
    m_timings.water    += msSince(t2, t3);
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
