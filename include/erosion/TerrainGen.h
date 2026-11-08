#pragma once

#include "erosion/AlignedAllocator.h"

#include <cstdint>

namespace erosion {

struct Grid;

/// Parametri fBm (fractional Brownian motion) generatora terena.
///
/// Koordinate šuma su normalizovane na [0, 1] preko mreže, tako da isti seed
/// daje isti oblik terena na 256², 512², 1024² i 2048² (samo s više detalja).
struct FbmSettings {
    uint32_t seed        = 42;
    float    baseFreq    = 4.0f; ///< Broj osnovnih "brda" preko cijele mreže
    int      octaves     = 6;
    float    persistence = 0.5f; ///< Faktor amplitude po oktavi
    float    lacunarity  = 2.0f; ///< Faktor frekvencije po oktavi
    float    exponent    = 1.0f; ///< h := h^exponent nakon normalizacije (>1 spljošti doline)
};

/// Gradient (Perlin-style) šum u tački (x, y). Rezultat je približno u [-1, 1].
/// Deterministički: zavisi samo od (x, y, seed), a ne od platforme.
[[nodiscard]] float gradientNoise(float x, float y, uint32_t seed);

/// fBm: zbir `octaves` slojeva gradientNoise, normalizovan na približno [-1, 1].
[[nodiscard]] float fbm(float x, float y, const FbmSettings& s);

/// Popunjava grid.b fBm terenom u opsegu [0, heightScale].
/// Ostala polja (voda, sediment, fluks, brzina) ostaju nepromijenjena.
void generateFbmTerrain(Grid& grid, const FbmSettings& settings, float heightScale);

/// Postavlja cijeli teren na konstantnu visinu (za testove: mirno jezero, dam-break).
void fillConstant(Field& field, float value);

/// Kopira normalizovanu heightmapu ([0, 1], npr. iz PGM-a) u grid.b, skaliranu na heightScale.
/// Dimenzije mape moraju odgovarati mreži.
void setTerrainFromNormalized(Grid& grid, const float* normalized, float heightScale);

} // namespace erosion
