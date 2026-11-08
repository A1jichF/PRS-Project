#pragma once

#include <cstddef>
#include <cstdint>
#include <new>
#include <vector>

namespace erosion {

namespace detail {

inline constexpr std::size_t kPageSize = 4096;
inline constexpr std::size_t kStagger  = 256; ///< 4 keš linije = 64 float-a
inline constexpr std::size_t kSlots    = kPageSize / kStagger; ///< 16 različitih pomaka

/// Sljedeći slot (pomak unutar stranice) za alokaciju u ovoj niti.
inline thread_local std::size_t tNextSlot = 0;

inline std::size_t nextStaggerOffset()
{
    return (tNextSlot++ % kSlots) * kStagger;
}

} // namespace detail

/// Određuje od kojeg slota počinju alokacije polja u ovom opsegu.
/// Grid::resize koristi slotove 0–8, Simulation za pomoćne bafere 9–11, pa
/// nijedna dva niza koja se koriste u istom kernelu nisu na istom pomaku.
class StaggerScope {
public:
    explicit StaggerScope(std::size_t firstSlot) noexcept
        : m_previous(detail::tNextSlot)
    {
        detail::tNextSlot = firstSlot;
    }
    ~StaggerScope() { detail::tNextSlot = m_previous; }

    StaggerScope(const StaggerScope&)            = delete;
    StaggerScope& operator=(const StaggerScope&) = delete;

private:
    std::size_t m_previous;
};

/// Alokator za polja mreže: poravnanje na `Align` bajta (64 B = keš linija,
/// pokriva i 32 B koje traži AVX2) i zadani pomak unutar 4 KB stranice.
///
/// Zašto pomak: polja iste mreže su iste veličine i alociraju se jedno za
/// drugim, pa bi inače bila na istom pomaku unutar stranice (npr. u i fR).
/// Kad petlja upiše u[i], a odmah zatim čita fR[i], procesor zbog poklapanja
/// nižih 12 bitova adrese pogrešno pretpostavi zavisnost i čeka ("4K aliasing").
/// Izmjereno na 256²: faza "voda" je bila 0.16, 0.38, 0.60 ili 0.82 ms/korak
/// zavisno od toga koliko se nizova slučajno poklopilo (README, "Mjerenje").
template <class T, std::size_t Align = 64>
struct AlignedAllocator {
    static_assert(Align <= detail::kStagger && detail::kStagger % Align == 0,
                  "pomak mora čuvati poravnanje");

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
        // Raspored: [ ...           | pokazivač na blok | podaci ]
        //           ^ blok poravnat na 4 KB              ^ blok + 4 KB + pomak
        const std::size_t offset = detail::kPageSize + detail::nextStaggerOffset();
        auto* base = static_cast<std::byte*>(
            ::operator new(n * sizeof(T) + offset, std::align_val_t{detail::kPageSize}));
        std::byte* data = base + offset;
        reinterpret_cast<std::byte**>(data)[-1] = base;
        return reinterpret_cast<T*>(data);
    }

    void deallocate(T* p, std::size_t) noexcept
    {
        std::byte* base = reinterpret_cast<std::byte**>(p)[-1];
        ::operator delete(base, std::align_val_t{detail::kPageSize});
    }

    template <class U>
    bool operator==(const AlignedAllocator<U, Align>&) const noexcept { return true; }
};

inline constexpr std::size_t kFieldAlignment = 64;

/// Jedno skalarno polje na mreži (npr. visina terena), poravnato za SIMD.
using Field = std::vector<float, AlignedAllocator<float, kFieldAlignment>>;

} // namespace erosion
