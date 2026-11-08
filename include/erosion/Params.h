#pragma once

namespace erosion {

/// Parametri simulacije (Mei, Decaudin, Hu 2007: virtual pipe model).
///
/// Konstante erozije (Kc, Ks, Kd, Ke) su početne vrijednosti i
/// podešavat ćemo ih u koracima 4–6, kada te faze budu implementirane.
struct Params {
    // ── Prostor i vrijeme ────────────────────────────────────────────────────
    /// Razmak ćelija u svjetskim jedinicama (= dužina virtualne cijevi l).
    /// 2.0 odgovara rendereru: 512 jedinica širine na 256 ćelija.
    float cellSize    = 2.0f;
    /// Maksimalna visina generisanog terena (kao heightScale u rendereru).
    float heightScale = 80.0f;
    float dt          = 0.02f;  ///< Vremenski korak [s]

    // ── Virtualne cijevi ─────────────────────────────────────────────────────
    float gravity  = 9.81f; ///< g
    float pipeArea = 1.0f;  ///< A, poprečni presjek cijevi

    // ── Voda ─────────────────────────────────────────────────────────────────
    float rainRate        = 0.0f;  ///< Ke_r: voda dodana po ćeliji u sekundi
    float evaporationRate = 0.015f;///< Ke: udio vode koji ispari u sekundi

    // ── Erozija i depozicija ─────────────────────────────────────────────────
    float sedimentCapacity = 1.0f;  ///< Kc
    float dissolvingRate   = 0.5f;  ///< Ks
    float depositionRate   = 1.0f;  ///< Kd
    float minTilt          = 0.05f; ///< Donja granica za sin(nagib), da C ne padne na 0 na ravnom

    // ── Termalna erozija (opcionalno) ────────────────────────────────────────
    float talusAngleDeg = 38.0f; ///< Ugao talusa
    float thermalRate   = 0.15f; ///< Udio viška materijala koji se prebaci po koraku
};

} // namespace erosion
