#ifndef EXACT_HUBBARD_SPECTRUM_HPP
#define EXACT_HUBBARD_SPECTRUM_HPP


#include <utility>
#include <vector>

#include "linalg.hpp"
#include "state.hpp"


/**
 * eigenstates are normalised
 * basis is / has to be normalised
 */
struct Spectrum
{
    IVector charges;
    DVector energies;
    std::vector<std::vector<std::size_t>> eigenStateIdxs;
    std::vector<std::vector<double>> eigenStateCoeffs;
    SumState basis;


    // ignores coeffs of basis
    static Spectrum compute(SumState const &inBasis);


    [[nodiscard]] std::size_t size() const noexcept;


private:
    explicit Spectrum(SumState const &inBasis);
};

#endif //EXACT_HUBBARD_SPECTRUM_HPP
