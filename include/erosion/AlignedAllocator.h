#pragma once

#include <cstddef>
#include <new>
#include <vector>

namespace erosion {

/// Alokator koji poravnava memoriju na `Align` bajta.
/// 64 B = jedna keš linija, a ujedno pokriva i 32 B koje traže AVX2 load/store.
template <class T, std::size_t Align = 64>
struct AlignedAllocator {
    using value_type = T;

    template <class U>
    struct rebind {
        using other = AlignedAllocator<U, Align>;
    };

    AlignedAllocator() noexcept = default;
    template <class U>
    AlignedAllocator(const AlignedAllocator<U, Align>&) noexcept {}

    [[nodiscard]] T* allocate(std::size_t n)
    {
        return static_cast<T*>(::operator new(n * sizeof(T), std::align_val_t{Align}));
    }

    void deallocate(T* p, std::size_t) noexcept
    {
        ::operator delete(p, std::align_val_t{Align});
    }

    template <class U>
    bool operator==(const AlignedAllocator<U, Align>&) const noexcept { return true; }
};

inline constexpr std::size_t kFieldAlignment = 64;

/// Jedno skalarno polje na mreži (npr. visina terena), poravnato za SIMD.
using Field = std::vector<float, AlignedAllocator<float, kFieldAlignment>>;

} // namespace erosion
