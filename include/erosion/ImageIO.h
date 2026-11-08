#pragma once

#include <string>
#include <vector>

namespace erosion {

/// Normalizovana slika visina: vrijednosti u [0, 1], red po red (y * width + x).
struct Heightmap {
    int width  = 0;
    int height = 0;
    std::vector<float> values;
};

/// Upisuje polje kao binarni PGM (P5). Vrijednosti se linearno mapiraju iz
/// [lo, hi] u [0, maxval] i odsijecaju izvan opsega.
/// bits16 = true daje 16-bitni PGM (maxval 65535), inače 8-bitni.
/// Vraća false ako se fajl ne može otvoriti.
bool writePGM(const std::string& path, const float* data, int width, int height,
              float lo, float hi, bool bits16 = true);

/// Isto kao gore, ali opseg [lo, hi] je min/max samih podataka.
bool writePGM(const std::string& path, const float* data, int width, int height,
              bool bits16 = true);

/// Čita binarni PGM (P5, 8 ili 16 bita) i vraća vrijednosti normalizovane na [0, 1].
/// Baca std::runtime_error ako fajl nije validan PGM.
[[nodiscard]] Heightmap readPGM(const std::string& path);

} // namespace erosion
