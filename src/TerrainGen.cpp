#include "erosion/TerrainGen.h"

#include "erosion/Grid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace erosion {

namespace {

/// Cjelobrojni hash (x, y, seed) -> 32 bita. Bez std::hash i bez <random>,
/// da rezultat bude isti na svim kompajlerima i platformama.
uint32_t hash2(int32_t x, int32_t y, uint32_t seed)
{
    uint32_t h = seed * 0x9E3779B9u;
    h ^= static_cast<uint32_t>(x) * 0x85EBCA6Bu;
    h ^= static_cast<uint32_t>(y) * 0xC2B2AE35u;
    h ^= h >> 16;
    h *= 0x7FEB352Du;
    h ^= h >> 15;
    h *= 0x846CA68Bu;
    h ^= h >> 16;
    return h;
}

// 8 jediničnih gradijenata (osi i dijagonale).
constexpr float kInvSqrt2 = 0.70710678f;
constexpr float kGrad[8][2] = {
    { 1.0f, 0.0f}, {-1.0f, 0.0f}, {0.0f, 1.0f}, {0.0f, -1.0f},
    { kInvSqrt2,  kInvSqrt2}, {-kInvSqrt2,  kInvSqrt2},
    { kInvSqrt2, -kInvSqrt2}, {-kInvSqrt2, -kInvSqrt2},
};

float gradDot(int32_t ix, int32_t iy, float dx, float dy, uint32_t seed)
{
    const float* g = kGrad[hash2(ix, iy, seed) & 7u];
    return g[0] * dx + g[1] * dy;
}

/// Quintic fade 6t^5 - 15t^4 + 10t^3 (Perlin 2002): glatki prvi i drugi izvod.
float fade(float t) { return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f); }

float lerp(float a, float b, float t) { return a + (b - a) * t; }

} // namespace

float gradientNoise(float x, float y, uint32_t seed)
{
    const float   fx = std::floor(x);
    const float   fy = std::floor(y);
    const int32_t ix = static_cast<int32_t>(fx);
    const int32_t iy = static_cast<int32_t>(fy);
    const float   dx = x - fx;
    const float   dy = y - fy;

    const float n00 = gradDot(ix,     iy,     dx,        dy,        seed);
    const float n10 = gradDot(ix + 1, iy,     dx - 1.0f, dy,        seed);
    const float n01 = gradDot(ix,     iy + 1, dx,        dy - 1.0f, seed);
    const float n11 = gradDot(ix + 1, iy + 1, dx - 1.0f, dy - 1.0f, seed);

    const float sx = fade(dx);
    const float sy = fade(dy);
    // 2D gradient šum je u [-sqrt(2)/2, sqrt(2)/2]; skaliramo na približno [-1, 1].
    return lerp(lerp(n00, n10, sx), lerp(n01, n11, sx), sy) * 1.41421356f;
}

float fbm(float x, float y, const FbmSettings& s)
{
    float sum      = 0.0f;
    float amp      = 1.0f;
    float freq     = 1.0f;
    float ampTotal = 0.0f;

    for (int o = 0; o < s.octaves; ++o) {
        // Različit seed po oktavi da se slojevi ne poklope na cjelobrojnim tačkama.
        sum      += amp * gradientNoise(x * freq, y * freq, s.seed + static_cast<uint32_t>(o) * 1013u);
        ampTotal += amp;
        amp      *= s.persistence;
        freq     *= s.lacunarity;
    }
    return ampTotal > 0.0f ? sum / ampTotal : 0.0f;
}

void generateFbmTerrain(Grid& grid, const FbmSettings& settings, float heightScale)
{
    if (grid.cellCount() == 0)
        throw std::invalid_argument("generateFbmTerrain: mreža nije alocirana");

    const int W = grid.width;
    const int H = grid.height;
    // Normalizujemo po većoj dimenziji, da ćelije ostanu kvadratne i u šumu.
    const float invExtent = 1.0f / static_cast<float>(std::max(W, H) - 1 > 0 ? std::max(W, H) - 1 : 1);

    float lo = std::numeric_limits<float>::max();
    float hi = std::numeric_limits<float>::lowest();

    for (int y = 0; y < H; ++y) {
        const float ny = static_cast<float>(y) * invExtent * settings.baseFreq;
        for (int x = 0; x < W; ++x) {
            const float nx = static_cast<float>(x) * invExtent * settings.baseFreq;
            const float h  = fbm(nx, ny, settings);
            grid.b[grid.idx(x, y)] = h;
            lo = std::min(lo, h);
            hi = std::max(hi, h);
        }
    }

    // Normalizacija na [0, 1], zatim oblikovanje (exponent) i skaliranje.
    const float range = hi - lo;
    const float inv   = range > 0.0f ? 1.0f / range : 0.0f;
    for (float& h : grid.b) {
        float t = (h - lo) * inv;
        if (settings.exponent != 1.0f)
            t = std::pow(t, settings.exponent);
        h = t * heightScale;
    }
}

void fillConstant(Field& field, float value)
{
    std::fill(field.begin(), field.end(), value);
}

void setTerrainFromNormalized(Grid& grid, const float* normalized, float heightScale)
{
    const std::size_t n = grid.cellCount();
    for (std::size_t i = 0; i < n; ++i)
        grid.b[i] = normalized[i] * heightScale;
}

} // namespace erosion
