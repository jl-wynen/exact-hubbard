#ifndef EXACT_HUBBARD_CONFIG_HPP
#define EXACT_HUBBARD_CONFIG_HPP

#include <array>
#include <cstdint>
#include <utility>


constexpr static std::size_t NSITES = 3;


using Link = std::pair<std::size_t, std::size_t>;

/*
 * If Link{i, j} is present, Link{j, i} is implicitly assumed to be there as well.
 */

// Triangle
constexpr static std::array nearestNeighbours
        = {Link{0, 1},
           Link{1, 2},
           Link{2, 0}};


constexpr double kappa = 1.0;

constexpr double beta = 1.0;

constexpr double U = 1.0;


#endif //EXACT_HUBBARD_CONFIG_HPP
