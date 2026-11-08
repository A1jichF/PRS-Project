#pragma once

namespace erosion {

/// Parametri simulacije (Mei, Decaudin, Hu 2007: virtual pipe model).
///
/// Razmak ćelija (dužina cijevi l) nije ovdje nego u Grid::cellSize, da ne
/// postoji na dva mjesta. Konstante erozije (Kc, Ks, Kd, Ke) su početne
/// vrijednosti i podešavat ćemo ih u koracima 4–6.
struct Params {
    float dt = 0.02f; ///< Vremenski korak [s]

    // ── Virtualne cijevi ─────────────────────────────────────────────────────
    float gravity  = 9.81f; ///< g
    /// A, poprečni presjek cijevi. Napomena: brzina talasa u pipe modelu ne
    /// zavisi od stvarne dubine vode, nego se ponaša kao da je dubina A / l.
    float pipeArea = 1.0f;

    // ── Voda ─────────────────────────────────────────────────────────────────
    float rainRate        = 0.0f;   ///< Visina vode dodana svakoj ćeliji u sekundi
    float evaporationRate = 0.015f; ///< Ke: udio vode koji ispari u sekundi (korak 6)
    /// Ispod ove prosječne dubine brzina se postavlja na 0, jer bi dijeljenje
    /// malom dubinom dalo ogromne (nefizikalne) brzine na skoro suhim ćelijama.
    float minDepthForVelocity = 1e-4f;

    // ── Erozija i depozicija (korak 4) ───────────────────────────────────────
    float sedimentCapacity = 1.0f;  ///< Kc
    float dissolvingRate   = 0.5f;  ///< Ks
    float depositionRate   = 1.0f;  ///< Kd
    float minTilt          = 0.05f; ///< Donja granica za sin(nagib), da C ne padne na 0 na ravnom

    // ── Termalna erozija (opcionalno) ────────────────────────────────────────
    float talusAngleDeg = 38.0f; ///< Ugao talusa
    float thermalRate   = 0.15f; ///< Udio viška materijala koji se prebaci po koraku
};

/// Parametri za pravljenje terena (ne utiču na samu simulaciju).
/// Zadane vrijednosti odgovaraju opengl-forward-rendereru: 512 jedinica
/// širine na 256 ćelija (cellSize = 2) i maksimalna visina 80.
struct WorldSettings {
    float cellSize    = 2.0f;
    float heightScale = 80.0f;
};

} // namespace erosion
