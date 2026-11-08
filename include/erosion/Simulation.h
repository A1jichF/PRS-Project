#pragma once

#include "erosion/Grid.h"
#include "erosion/Kernels.h"
#include "erosion/Params.h"

#include <vector>

namespace erosion {

/// Vrijeme po fazi, sabrano preko svih koraka od zadnjeg resetTimings() [ms].
struct PhaseTimings {
    double addWater = 0.0;
    double flux     = 0.0;
    double water    = 0.0;
    double erosion  = 0.0;
    long   steps    = 0;

    [[nodiscard]] double total() const { return addWater + flux + water + erosion; }
};

/// Vodi simulaciju: drži stanje mreže, parametre i izvore vode, i redom
/// poziva kernele za svaki korak.
///
/// Redoslijed faza u jednom koraku (Mei et al. 2007):
///   1. dodavanje vode              ✔ korak 2
///   2. fluks kroz virtualne cijevi ✔ korak 2
///   3. visina vode + brzina        ✔ korak 2
///   4. erozija i depozicija        ✔ korak 4
///   5. transport sedimenta           korak 5
///   6. isparavanje                   korak 6
class Simulation {
public:
    explicit Simulation(Grid grid, Params params = {});

    /// Jedan vremenski korak (sve implementirane faze).
    void step();

    /// n koraka zaredom.
    void run(int n);

    void addSource(const WaterSource& source) { m_sources.push_back(source); }
    void clearSources() { m_sources.clear(); }

    [[nodiscard]] Grid&                     grid() { return m_grid; }
    [[nodiscard]] const Grid&               grid() const { return m_grid; }
    [[nodiscard]] Params&                   params() { return m_params; }
    [[nodiscard]] const Params&             params() const { return m_params; }
    [[nodiscard]] const PhaseTimings&       timings() const { return m_timings; }
    [[nodiscard]] long                      stepCount() const { return m_stepCount; }
    [[nodiscard]] double                    time() const { return m_time; }

    void resetTimings() { m_timings = {}; }

private:
    Grid                     m_grid;
    Params                   m_params;
    Field                    m_bScratch; ///< Drugi bafer terena za erodeDeposit
    std::vector<WaterSource> m_sources;
    PhaseTimings             m_timings;
    long                     m_stepCount = 0;
    double                   m_time      = 0.0; ///< Simulirano vrijeme [s]
};

} // namespace erosion
