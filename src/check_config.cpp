#include "config.hpp"


/**
 * Compiling this file checks whether settings in config.hpp are consistent.
 */

static_assert(NSITES > 0, "There must be more than 0 sites.");


namespace {
    template <std::size_t N>
    constexpr bool sitesAreInRange(std::array<Link, N> const &links) noexcept
    {
        for (auto const [a, b] : links) {
            if (a >= NSITES or b >= NSITES) {
                return false;
            }
        }
        return true;
    }
}

static_assert(sitesAreInRange(nearestNeighbours),
              "All sites in nearestNeighbours must be between 0 and NSITES");


namespace {
    template <std::size_t N>
    constexpr bool containsEverySite(std::array<Link, N> const &links) noexcept
    {
        std::array<std::size_t, NSITES> counts{};
        for (auto const link : links) {
            counts[link.first]++;
            counts[link.second]++;
        }

        for (std::size_t i = 0; i < NSITES; ++i) {
            if (counts[i] == 0) {
                return false;
            }
        }

        return true;
    }
}

static_assert(containsEverySite(nearestNeighbours),
              "nearestNeighbours must contain links for every site on the lattice.");
