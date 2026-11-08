#pragma once

// Kerneli simulacije, jedan po fazi algoritma (Mei et al. 2007, poglavlje 3).
//
// Namespace `seq` je sekvencijalna (referentna) verzija. OpenMP i SIMD verzije
// dobit će svoje namespace-ove s istim potpisima, tako da ih Simulation može
// mijenjati i da testovi mogu porediti rezultate.

#include "erosion/AlignedAllocator.h"

#include <span>

namespace erosion {

struct Grid;
struct Params;

/// Izvor vode: dodaje `rate` visine vode u sekundi svakoj ćeliji unutar kruga.
struct WaterSource {
    int   x      = 0;
    int   y      = 0;
    float radius = 1.0f;
    float rate   = 1.0f;
};

namespace seq {

/// Faza 1: dodavanje vode (ravnomjerna kiša Params::rainRate + izvori).
///   d1 = d + dt * r
void addWater(Grid& g, const Params& p, std::span<const WaterSource> sources);

/// Faza 2: odlazni fluks kroz virtualne cijevi prema 4 susjeda.
///   f_new = max(0, f + dt * A * g * Δh / l),  Δh = (b + d)_ćelija - (b + d)_susjed
/// Zatim se sva 4 fluksa skaliraju s K = min(1, d * l² / (Σf * dt)), da ćelija
/// ne pošalje više vode nego što ima. Rubovi mreže su zatvoreni zidovi (fluks 0).
///
/// Radi "u mjestu": čita b i d susjeda (koji se ovdje ne mijenjaju) i svoj stari fluks.
void computeFlux(Grid& g, const Params& p);

/// Faza 3: nova visina vode i polje brzine.
///   d2 = d1 + dt * (Σ dolazni fluks - Σ odlazni fluks) / l²
///   u  = ΔWX / (l * d̄),  ΔWX = (fR_lijevog - fL + fR - fL_desnog) / 2,  d̄ = (d1 + d2) / 2
///   (analogno v po y osi)
///
/// Radi "u mjestu": čita fluks susjeda (ovdje se ne mijenja) i svoje d.
void updateWaterAndVelocity(Grid& g, const Params& p);

/// Faza 4: erozija i depozicija.
///   sin α  iz centralnih razlika terena (na rubu jednostrano)
///   C      = Kc * max(sin α, minTilt) * |(u, v)| * lmax(d)
///   C > s: b -= Ks*dt*(C - s), s += isto     (erozija)
///   C ≤ s: b += Kd*dt*(s - C), s -= isto     (depozicija)
///
/// Double-buffering: nagib se računa iz STAROG terena susjeda, pa se novi teren
/// piše u `bNext`, a na kraju se g.b i bNext zamijene (O(1), samo pokazivači).
/// Tako rezultat ne zavisi od redoslijeda obrade ćelija. `bNext` mora imati
/// g.cellCount() elemenata; nakon poziva sadrži stari teren.
/// Sediment s se mijenja u mjestu (čita se samo vlastita ćelija).
void erodeDeposit(Grid& g, const Params& p, Field& bNext);

} // namespace seq
} // namespace erosion
