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

// Triangle
[[maybe_unused]] constexpr static std::array nearestNeighbours
        = {Link{0, 1},
           Link{1, 2},
           Link{2, 0}};

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
constexpr double U = 1.0;

/// Number of lattice sites.
constexpr static std::size_t NSITES = nearestNeighbours.size();


#endif //EXACT_HUBBARD_CONFIG_HPP
