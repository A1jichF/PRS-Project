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
/// Stijena (Grid::r): erozija uzima najviše sloj tla b - r. Ostatak tražene
/// količine se uzima iz stijene pomnožen s Params::rockErodibility (0 = ništa),
/// i tada se spuštaju i b i r. Depozicija ne dira r: nataloženi materijal je
/// novo tlo iznad stijene. Uvijek vrijedi b >= r. Polje r se mijenja u mjestu
/// (čita se i piše samo vlastita ćelija).
///
/// Double-buffering: nagib se računa iz STAROG terena susjeda, pa se novi teren
/// piše u `bNext`, a na kraju se g.b i bNext zamijene (O(1), samo pokazivači).
/// Tako rezultat ne zavisi od redoslijeda obrade ćelija. `bNext` mora imati
/// g.cellCount() elemenata; nakon poziva sadrži stari teren.
/// Sediment s se mijenja u mjestu (čita se samo vlastita ćelija).
void erodeDeposit(Grid& g, const Params& p, Field& bNext);

/// Faza 5: transport sedimenta, semi-Lagrangian advekcija.
///   s_new(x, y) = s(x - u*dt/l, y - v*dt/l)   (u ćelijama, bilinearna interpolacija)
/// Tačka iz koje se uzima sediment odsijeca se na granice mreže.
///
/// Double-buffering: čita se STARI s (i susjeda), novi se piše u `sNext`, pa se
/// g.s i sNext zamijene. Pristup memoriji zavisi od brzine (nepravilan "gather"),
/// što je važno za kasniju SIMD/OpenMP analizu.
///
/// Napomena: metoda NE čuva masu. Ona rješava ∂s/∂t + v·∇s = 0, a za plitku
/// vodu bi trebalo ∂s/∂t + ∇·(s v) = 0. Razlika je s·(∇·v): gdje se tok
/// skuplja (∇·v < 0, npr. ulaz u jezero), sediment nestaje.
void advectSediment(Grid& g, const Params& p, Field& sNext);

/// Faza 5 (alternativa, zadana): konzervativni transport sedimenta kroz cijevi.
///
/// Poziva se POSLIJE computeFlux, a PRIJE updateWaterAndVelocity, dok g.d još
/// sadrži d1 (vodu prije toka). Ćelija i šalje susjedu j udio svog sedimenta
/// jednak udjelu vode koji šalje:
///   r_i     = dt / max(d1_i * l², Σf_i * dt)      (udio po jedinici fluksa, Σ udjela <= 1)
///   s_new_i = s_i * (1 - Σf_i * r_i) + Σ_j s_j * f_{j→i} * r_j
/// Isti r_j se koristi i kad j računa odlazni i kad i računa dolazni sediment,
/// pa se ukupna masa čuva do greške zaokruživanja. Ako je s = c*d svuda
/// (ista koncentracija), poslije koraka vrijedi s = c*d2: sediment prati vodu.
///
/// `sNext` i `rate` su pomoćni baferi veličine g.cellCount(); poslije poziva
/// sNext sadrži stari sediment (zamjena bafera).
void transportSedimentFlux(Grid& g, const Params& p, Field& sNext, Field& rate);

/// Faza 6: isparavanje.
///   d = d * max(0, 1 - Ke * dt)
/// Isparava samo voda; sediment ostaje i, kako voda nestaje, taloži se kroz
/// fazu 4 (kapacitet C pada s brzinom i dubinom).
void evaporate(Grid& g, const Params& p);

} // namespace seq
} // namespace erosion
