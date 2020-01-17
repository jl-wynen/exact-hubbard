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



/**
 * Turn matrix elements for an operator in basis `spectrum.basis`
 * into matrix elements in the basis of eigenvectors.
 * \param matrix Matrix elements in `spectrum.basis`.
 * \param spectrum Provides basis and eigenstates.
 * \return Matrix elements in eigenbasis.
 */
DSparseMatrix toEigenspaceMatrix(DMatrix const &matrix, Spectrum const &spectrum);


#endif //EXACT_HUBBARD_SPECTRUM_HPP
