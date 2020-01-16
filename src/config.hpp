#ifndef EXACT_HUBBARD_CONFIG_HPP
#define EXACT_HUBBARD_CONFIG_HPP

/**
 * Configure the program.
 * Specify parameters and lattice geometry in this file and compile.
 */


#include <array>
#include <cstdint>
#include <utility>


/// Type for links in nearest-neighbour graph.
using Link = std::pair<std::size_t, std::size_t>;

/*
 * Encode Lattice geometry by specifying nearest-neighbour relations.
 *
 * The lattice is assumed to be symmetric.
 * If Link{i, j} is specified, Link{j, i} is implicitly assumed
 * to be there as well.
 */

// Two sites
//[[maybe_unused]] constexpr static std::array nearestNeighbours
//        = {Link{0, 1}};

// Triangle
[[maybe_unused]] constexpr static std::array nearestNeighbours
        = {Link{0, 1},
           Link{1, 2},
           Link{2, 0}};

// Square, needs kappa=2
//[[maybe_unused]] constexpr static std::array nearestNeighbours
//        = {Link{0, 1},
//           Link{0, 3},
//           Link{1, 2},
//           Link{2, 3}};

// Pentagon
//[[maybe_unused]] constexpr static std::array nearestNeighbours
//        = {Link{0, 1},
//           Link{1, 2},
//           Link{2, 3},
//           Link{3, 4},
//           Link{4, 0}};


/// Nearest-neighbour hopping parameter.
constexpr double kappa = 1.0;

/// Inverse temperature.
constexpr double beta = 1.0;

/// On-site interaction strength.
constexpr double U = 4.0;

/// Number of time slices.
constexpr std::size_t NT = 32;


/// Compute the number of lattice sites from nearestNeighbours.
constexpr std::size_t computeNumSites()
{
    std::size_t nsites = 0;
    for (auto const [a, b] : nearestNeighbours) {
        if (a+1 > nsites) {
            nsites = a + 1;
        }
        if (b+1 > nsites) {
            nsites = b + 1;
        }
    }
    return nsites;
}


/// Number of lattice sites.
constexpr static std::size_t NSITES = computeNumSites();


#endif //EXACT_HUBBARD_CONFIG_HPP
