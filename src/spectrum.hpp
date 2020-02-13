#ifndef EXACT_HUBBARD_SPECTRUM_HPP
#define EXACT_HUBBARD_SPECTRUM_HPP

/** \file
 * \brief Spectrum storage and computation.
 */

#include <utility>
#include <vector>

#include "linalg.hpp"
#include "state.hpp"


/**
 * Stores an energy spectrum and associated eigenstates.
 *
 * The eigenstates are simultaneous eigenvectors of the Hamiltonian and charge operators.
 * The basis and all eigenstates are normalised.
 *
 * The states are stored via index and coefficient lists based on the basis.
 * You can obtain eigenstate `i` via
 * ```{.cpp}
// construct eigenstate i
SumState state;
auto const &indexes = spectrum.eigenStateIdxs[i];
auto const &coeffs = spectrum.eigenStateCoeffs[i];
for (std::size_t j = 0; j < indexes.size(); ++j) {
    auto const [c, e] = spectrum.basis[indexes[j]];
    state.push(coeffs[j]*c, e);
}
 * ```
 */
struct Spectrum
{
    /// Expectation value of the charge operator for each eigenstate.
    IVector charges;
    /// Expectation value of the Hamiltonian for each eigenstate.
    DVector energies;
    /// Indexes of eigenstates into the basis.
    std::vector<std::vector<std::size_t>> eigenStateIdxs;
    /// Coefficients of eigenstates.
    std::vector<std::vector<double>> eigenStateCoeffs;
    /// Basis elements.
    SumState basis;


    /// Computes the spectrum for a given basis.
    /**
     * \param inBasis Basis states, must be normalised.
     * \return A new instance of Spectrum-
     */
    static Spectrum compute(SumState const &inBasis);


    /// Return the number of eigenstates.
    [[nodiscard]] std::size_t size() const noexcept;


private:
    /// Disallow construction from the outside.
    explicit Spectrum(SumState const &inBasis);
};



/**
 * Turn matrix elements of an operator in basis `spectrum.basis`
 * into matrix elements in the basis of eigenvectors.
 * \param matrix Matrix elements in `spectrum.basis`.
 * \param spectrum Provides basis and eigenstates.
 * \return Matrix elements in eigenbasis.
 */
DSparseMatrix toEigenspaceMatrix(DMatrix const &matrix, Spectrum const &spectrum);


#endif //EXACT_HUBBARD_SPECTRUM_HPP
