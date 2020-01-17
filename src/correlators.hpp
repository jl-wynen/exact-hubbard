#ifndef EXACT_HUBBARD_CORRELATORS_HPP
#define EXACT_HUBBARD_CORRELATORS_HPP

/** \file
 * \brief Correlator storage and computation.
 */

#include <cassert>
#include <vector>

#include "config.hpp"
#include "spectrum.hpp"


struct Correlators
{
    std::vector<double> data;


    explicit Correlators()
            : data(NSITES * NSITES * NT)
    { }


    [[nodiscard]] static std::size_t totalIndex(std::size_t const i,
                                                std::size_t const j,
                                                std::size_t const t) noexcept
    {
        assert(i < NSITES);
        assert(j < NSITES);
        assert(t < NT);
        return (i*NSITES + j)*NT + t;
    }


    double operator()(std::size_t const i,
                      std::size_t const j,
                      std::size_t const t) const noexcept
    {
        return data[totalIndex(i, j, t)];
    }


    double &operator()(std::size_t const i,
                       std::size_t const j,
                       std::size_t const t) noexcept
    {
        return data[totalIndex(i, j, t)];
    }

};


Correlators computeCorrelators(Spectrum const &spectrum);

#endif //EXACT_HUBBARD_CORRELATORS_HPP
