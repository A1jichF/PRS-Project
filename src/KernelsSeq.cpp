#include "erosion/Grid.h"
#include "erosion/Kernels.h"
#include "erosion/Params.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace erosion::seq {

void addWater(Grid& g, const Params& p, std::span<const WaterSource> sources)
{
    if (p.rainRate != 0.0f) {
        const float add = p.dt * p.rainRate;
        for (float& d : g.d)
            d += add;
    }

    for (const WaterSource& src : sources) {
        const float add = p.dt * src.rate;
        const float r2  = src.radius * src.radius;
        const int   r   = static_cast<int>(src.radius) + 1;
        const int   x0  = std::max(0, src.x - r);
        const int   x1  = std::min(g.width - 1, src.x + r);
        const int   y0  = std::max(0, src.y - r);
        const int   y1  = std::min(g.height - 1, src.y + r);
        for (int y = y0; y <= y1; ++y) {
            for (int x = x0; x <= x1; ++x) {
                const float dx = static_cast<float>(x - src.x);
                const float dy = static_cast<float>(y - src.y);
                if (dx * dx + dy * dy <= r2)
                    g.d[g.idx(x, y)] += add;
            }
        }
    }
}

void computeFlux(Grid& g, const Params& p)
{
    const int    W        = g.width;
    const int    H        = g.height;
    const float  l        = g.cellSize;
    const float  cellArea = l * l;
    const float  k        = p.dt * p.pipeArea * p.gravity / l;
    const float* b        = g.b.data();
    const float* d        = g.d.data();
    float*       fL       = g.fL.data();
    float*       fR       = g.fR.data();
    float*       fT       = g.fT.data();
    float*       fB       = g.fB.data();

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            const std::size_t i = g.idx(x, y);
            const float       h = b[i] + d[i];

            // Zatvoreni rubovi: prema van mreže nema cijevi.
            float l_ = 0.0f, r_ = 0.0f, t_ = 0.0f, b_ = 0.0f;
            if (x > 0)     l_ = std::max(0.0f, fL[i] + k * (h - (b[i - 1] + d[i - 1])));
            if (x < W - 1) r_ = std::max(0.0f, fR[i] + k * (h - (b[i + 1] + d[i + 1])));
            if (y > 0)     t_ = std::max(0.0f, fT[i] + k * (h - (b[i - W] + d[i - W])));
            if (y < H - 1) b_ = std::max(0.0f, fB[i] + k * (h - (b[i + W] + d[i + W])));

            // Ćelija ne smije poslati više vode nego što ima: Σf * dt <= d * l².
            const float sum = l_ + r_ + t_ + b_;
            if (sum > 0.0f) {
                const float K = std::min(1.0f, d[i] * cellArea / (sum * p.dt));
                l_ *= K;
                r_ *= K;
                t_ *= K;
                b_ *= K;
            }

            fL[i] = l_;
            fR[i] = r_;
            fT[i] = t_;
            fB[i] = b_;
        }
    }
}

void updateWaterAndVelocity(Grid& g, const Params& p)
{
    const int    W           = g.width;
    const int    H           = g.height;
    const float  l           = g.cellSize;
    const float  invCellArea = 1.0f / (l * l);
    const float* fL          = g.fL.data();
    const float* fR          = g.fR.data();
    const float* fT          = g.fT.data();
    const float* fB          = g.fB.data();
    float*       d           = g.d.data();
    float*       u           = g.u.data();
    float*       v           = g.v.data();

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            const std::size_t i = g.idx(x, y);

            // Dolazni fluks: šta susjedi šalju prema ovoj ćeliji.
            const float inL = x > 0     ? fR[i - 1] : 0.0f; // lijevi susjed šalje desno
            const float inR = x < W - 1 ? fL[i + 1] : 0.0f; // desni susjed šalje lijevo
            const float inT = y > 0     ? fB[i - W] : 0.0f; // gornji susjed šalje dolje
            const float inB = y < H - 1 ? fT[i + W] : 0.0f; // donji susjed šalje gore

            const float inflow  = inL + inR + inT + inB;
            const float outflow = fL[i] + fR[i] + fT[i] + fB[i];

            const float d1 = d[i];
            // K-skaliranje garantuje d2 >= 0 u egzaktnoj aritmetici; max uklanja
            // samo grešku zaokruživanja reda 1e-7.
            const float d2 = std::max(0.0f, d1 + p.dt * (inflow - outflow) * invCellArea);
            d[i] = d2;

            // Prosječan protok kroz ćeliju po x i y (Mei et al., jednačine 8–9).
            const float dWX  = 0.5f * (inL - fL[i] + fR[i] - inR);
            const float dWY  = 0.5f * (inT - fT[i] + fB[i] - inB);
            const float dBar = 0.5f * (d1 + d2);

            if (dBar > p.minDepthForVelocity) {
                const float inv = 1.0f / (l * dBar);
                u[i] = dWX * inv;
                v[i] = dWY * inv;
            } else {
                u[i] = 0.0f;
                v[i] = 0.0f;
            }
        }
    }
}

void erodeDeposit(Grid& g, const Params& p, Field& bNext)
{
    const int    W  = g.width;
    const int    H  = g.height;
    const float  l  = g.cellSize;
    const float  ks = std::min(1.0f, p.dissolvingRate * p.dt);
    const float  kd = std::min(1.0f, p.depositionRate * p.dt);
    const float  invMaxDepth = p.maxErosionDepth > 0.0f ? 1.0f / p.maxErosionDepth : 0.0f;
    const float* b  = g.b.data();
    const float* d  = g.d.data();
    const float* u  = g.u.data();
    const float* v  = g.v.data();
    float*       s  = g.s.data();
    float*       bn = bNext.data();

    for (int y = 0; y < H; ++y) {
        // Centralna razlika u unutrašnjosti, jednostrana na rubu.
        const int   yu   = y > 0 ? y - 1 : y;
        const int   yd   = y < H - 1 ? y + 1 : y;
        const float invY = yd > yu ? 1.0f / (static_cast<float>(yd - yu) * l) : 0.0f;

        for (int x = 0; x < W; ++x) {
            const int   xl   = x > 0 ? x - 1 : x;
            const int   xr   = x < W - 1 ? x + 1 : x;
            const float invX = xr > xl ? 1.0f / (static_cast<float>(xr - xl) * l) : 0.0f;

            const std::size_t i  = g.idx(x, y);
            const float       gx = (b[g.idx(xr, y)] - b[g.idx(xl, y)]) * invX;
            const float       gy = (b[g.idx(x, yd)] - b[g.idx(x, yu)]) * invY;

            // sin α = |∇b| / sqrt(1 + |∇b|²)
            const float grad2 = gx * gx + gy * gy;
            const float sinA  = std::max(p.minTilt, std::sqrt(grad2 / (1.0f + grad2)));
            const float speed = std::sqrt(u[i] * u[i] + v[i] * v[i]);
            const float lmax  = invMaxDepth > 0.0f ? std::min(1.0f, d[i] * invMaxDepth) : 1.0f;
            const float C     = p.sedimentCapacity * sinA * speed * lmax;

            if (C > s[i]) {
                const float amount = ks * (C - s[i]);
                bn[i] = b[i] - amount;
                s[i] += amount;
            } else {
                const float amount = kd * (s[i] - C);
                bn[i] = b[i] + amount;
                s[i] -= amount;
            }
        }
    }

    std::swap(g.b, bNext);
}

void advectSediment(Grid& g, const Params& p, Field& sNext)
{
    const int    W     = g.width;
    const int    H     = g.height;
    const float  scale = p.dt / g.cellSize; // brzina [jedinica/s] -> pomak u ćelijama
    const float  maxX  = static_cast<float>(W - 1);
    const float  maxY  = static_cast<float>(H - 1);
    const float* s     = g.s.data();
    const float* u     = g.u.data();
    const float* v     = g.v.data();
    float*       sn    = sNext.data();

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            const std::size_t i = g.idx(x, y);

            // Odakle je voda došla u ovu ćeliju (unazad po brzini).
            const float fx = std::clamp(static_cast<float>(x) - u[i] * scale, 0.0f, maxX);
            const float fy = std::clamp(static_cast<float>(y) - v[i] * scale, 0.0f, maxY);

            const int   x0 = static_cast<int>(fx); // fx >= 0, pa je ovo floor
            const int   y0 = static_cast<int>(fy);
            const int   x1 = std::min(x0 + 1, W - 1);
            const int   y1 = std::min(y0 + 1, H - 1);
            const float tx = fx - static_cast<float>(x0);
            const float ty = fy - static_cast<float>(y0);

            const float s00 = s[g.idx(x0, y0)];
            const float s10 = s[g.idx(x1, y0)];
            const float s01 = s[g.idx(x0, y1)];
            const float s11 = s[g.idx(x1, y1)];

            const float top    = s00 + (s10 - s00) * tx;
            const float bottom = s01 + (s11 - s01) * tx;
            sn[i] = top + (bottom - top) * ty;
        }
    }

    std::swap(g.s, sNext);
}

void transportSedimentFlux(Grid& g, const Params& p, Field& sNext, Field& rate)
{
    const int    W        = g.width;
    const int    H        = g.height;
    const float  cellArea = g.cellSize * g.cellSize;
    const float  dt       = p.dt;
    const float* d        = g.d.data(); // d1: voda prije ažuriranja
    const float* s        = g.s.data();
    const float* fL       = g.fL.data();
    const float* fR       = g.fR.data();
    const float* fT       = g.fT.data();
    const float* fB       = g.fB.data();
    float*       r        = rate.data();
    float*       sn       = sNext.data();
    const std::size_t n   = g.cellCount();

    // Ispod ovog odlaznog volumena ćelija ne šalje sediment. Bez praga bi za
    // fluks reda 1e-40 (denormalni broj) r = dt / (out*dt) prešao opseg float-a
    // (inf), a 0 * inf = NaN. Obje strane koriste isti r, pa se masa i dalje čuva.
    constexpr float kMinOutVolume = 1e-20f;

    // Prolaz 1: udio vode (pa i sedimenta) po jedinici fluksa za svaku ćeliju.
    for (std::size_t i = 0; i < n; ++i) {
        const float out    = fL[i] + fR[i] + fT[i] + fB[i];
        const float outVol = out * dt;
        r[i] = outVol > kMinOutVolume ? dt / std::max(d[i] * cellArea, outVol) : 0.0f;
    }

    // Prolaz 2: ono što ostane + ono što susjedi pošalju (gather, bez konflikta pisanja).
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            const std::size_t i   = g.idx(x, y);
            const float       out = fL[i] + fR[i] + fT[i] + fB[i];

            // out * r <= 1 u egzaktnoj aritmetici; zaokruživanje može dati
            // 1 + 1e-7, pa bi ćelija koja šalje svu vodu dobila negativan sediment.
            float value = s[i] * std::max(0.0f, 1.0f - out * r[i]);
            if (x > 0)     value += s[i - 1] * fR[i - 1] * r[i - 1]; // lijevi šalje desno
            if (x < W - 1) value += s[i + 1] * fL[i + 1] * r[i + 1]; // desni šalje lijevo
            if (y > 0)     value += s[i - W] * fB[i - W] * r[i - W]; // gornji šalje dolje
            if (y < H - 1) value += s[i + W] * fT[i + W] * r[i + W]; // donji šalje gore
            sn[i] = value;
        }
    }

    std::swap(g.s, sNext);
}

void evaporate(Grid& g, const Params& p)
{
    if (p.evaporationRate == 0.0f)
        return;
    // max(0, ...) da prevelik Ke * dt ne napravi negativnu vodu.
    const float factor = std::max(0.0f, 1.0f - p.evaporationRate * p.dt);
    for (float& d : g.d)
        d *= factor;
}

} // namespace erosion::seq
