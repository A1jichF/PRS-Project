#include "erosion/ImageIO.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace erosion {

bool writePGM(const std::string& path, const float* data, int width, int height,
              float lo, float hi, bool bits16)
{
    std::ofstream out(path, std::ios::binary);
    if (!out)
        return false;

    const int maxval = bits16 ? 65535 : 255;
    out << "P5\n" << width << ' ' << height << '\n' << maxval << '\n';

    const float range = hi - lo;
    const float scale = range > 0.0f ? static_cast<float>(maxval) / range : 0.0f;
    const std::size_t n = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);

    std::vector<unsigned char> bytes(n * (bits16 ? 2 : 1));
    for (std::size_t i = 0; i < n; ++i) {
        float q = (data[i] - lo) * scale;
        q = std::clamp(q, 0.0f, static_cast<float>(maxval));
        const auto value = static_cast<uint32_t>(std::lround(q));
        if (bits16) {
            // PGM standard: 16-bitne vrijednosti su big-endian.
            bytes[2 * i]     = static_cast<unsigned char>(value >> 8);
            bytes[2 * i + 1] = static_cast<unsigned char>(value & 0xFFu);
        } else {
            bytes[i] = static_cast<unsigned char>(value);
        }
    }
    out.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    return static_cast<bool>(out);
}

bool writePGM(const std::string& path, const float* data, int width, int height, bool bits16)
{
    const std::size_t n = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    float lo = 0.0f;
    float hi = 0.0f;
    if (n > 0) {
        const auto [mn, mx] = std::minmax_element(data, data + n);
        lo = *mn;
        hi = *mx;
    }
    return writePGM(path, data, width, height, lo, hi, bits16);
}

namespace {

/// Čita sljedeći broj iz PGM zaglavlja, preskačući razmake i komentare (#...).
int readHeaderInt(std::istream& in)
{
    for (;;) {
        const int c = in.peek();
        if (c == '#') {
            in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } else if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            in.get();
        } else {
            break;
        }
    }
    int value = 0;
    if (!(in >> value))
        throw std::runtime_error("readPGM: neispravno zaglavlje");
    return value;
}

} // namespace

Heightmap readPGM(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
        throw std::runtime_error("readPGM: ne mogu otvoriti " + path);

    char magic[2] = {};
    in.read(magic, 2);
    if (!in || magic[0] != 'P' || magic[1] != '5')
        throw std::runtime_error("readPGM: podržan je samo binarni PGM (P5): " + path);

    Heightmap hm;
    hm.width  = readHeaderInt(in);
    hm.height = readHeaderInt(in);
    const int maxval = readHeaderInt(in);
    if (hm.width <= 0 || hm.height <= 0 || maxval <= 0 || maxval > 65535)
        throw std::runtime_error("readPGM: neispravne dimenzije ili maxval: " + path);
    in.get(); // tačno jedan razmak između zaglavlja i podataka

    const bool        wide = maxval > 255;
    const std::size_t n    = static_cast<std::size_t>(hm.width) * static_cast<std::size_t>(hm.height);
    std::vector<unsigned char> bytes(n * (wide ? 2 : 1));
    in.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    if (in.gcount() != static_cast<std::streamsize>(bytes.size()))
        throw std::runtime_error("readPGM: fajl je kraći od očekivanog: " + path);

    hm.values.resize(n);
    const float inv = 1.0f / static_cast<float>(maxval);
    for (std::size_t i = 0; i < n; ++i) {
        const uint32_t raw = wide ? (static_cast<uint32_t>(bytes[2 * i]) << 8) | bytes[2 * i + 1]
                                  : bytes[i];
        hm.values[i] = static_cast<float>(raw) * inv;
    }
    return hm;
}

} // namespace erosion
