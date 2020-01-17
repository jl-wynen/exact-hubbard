#include "spectrum.hpp"

#include <blaze/math/lapack/syev.h>
#include <blaze/math/Submatrix.h>

#include "operator.hpp"


namespace {
    class EqualChargeIter
    {
        SumState basis_;
        std::size_t currentI_ = 0;
        int currentCharge_ = 0;
        ChargeOperator Q_{};


    public:
        explicit EqualChargeIter(SumState basis) : basis_{std::move(basis)} { }


        [[nodiscard]] bool finished() const noexcept
        {
            return currentI_ == basis_.size();
        }


        std::pair<SumState, int> next()
        {
            if (currentI_ == basis_.size()) {
                throw std::runtime_error("EqualChargeIter is out of states");
            }

            SumState res;
            currentCharge_ = Q_.computeCharge(basis_[currentI_].second);
            res.push(basis_[currentI_].first, basis_[currentI_].second);
            currentI_++;

            for (; currentI_ < basis_.size()
                   and Q_.computeCharge(basis_[currentI_].second) == currentCharge_;
                   ++currentI_)
            {
                res.push(basis_[currentI_].first, basis_[currentI_].second);
            }

            return std::pair{res, currentCharge_};
        };
    };


    template <typename Vec>
    SumState stateInBasis(Vec const &coefs, SumState const &basis)
    {
        SumState estate;
        for (std::size_t i = 0; i < coefs.size(); ++i) {
            if (std::abs(coefs[i]) > 1e-13) {
                auto const [basisCoef, basisState] = basis[i];
                estate.push(coefs[i] * basisCoef, basisState);
            }
        }
        return estate;
    }


    void computeSubSpectrum(SumState const &basis, int const charge, Spectrum &out, std::size_t &insertionOffset)
    {
        // compute spectrum
        SumOperator hamiltonian{ParticleHop{},
                                HoleHop{},
                                SquaredNumberOperator{}};
        DMatrix matrix = toMatrix(hamiltonian, basis);
        DVector evals(matrix.rows());
        blaze::syev(matrix, evals, 'V', 'U');

        // store spectrum
        for (std::size_t i = 0; i < evals.size(); ++i) {
            out.charges[insertionOffset + i] = charge;
            out.energies[insertionOffset + i] = evals[i];

            // `syev` stores the eigenvectors row-wise in `matrix`.
            for (std::size_t j = 0; j < matrix.columns(); ++j) {
                if (double const coef = matrix(i, j); std::abs(coef) > 1e-13) {
                    out.eigenStateIdxs[insertionOffset + i].push_back(insertionOffset + j);
                    out.eigenStateCoeffs[insertionOffset + i].push_back(coef);
                }
            }
        }

        insertionOffset += evals.size();
    }
}


Spectrum::Spectrum(SumState const &inBasis)
        : charges(inBasis.size()), energies(inBasis.size()),
          eigenStateIdxs(inBasis.size()), eigenStateCoeffs(inBasis.size()),
          basis(inBasis)
{ }


Spectrum Spectrum::compute(SumState const &inBasis)
{
    Spectrum spectrum(inBasis);

    // sort wrt. charge
    ChargeOperator Q{};
    std::sort(spectrum.basis.states(), spectrum.basis.states() + spectrum.basis.size(),
              [&Q](State const &a, State const &b) {
                  return Q.computeCharge(a) < Q.computeCharge(b);
              });

    // compute spectrum for given charge
    std::size_t insertionOffset = 0;
    for (EqualChargeIter eci{spectrum.basis}; not eci.finished();) {
        auto const [subBasis, charge] = eci.next();
        computeSubSpectrum(subBasis, charge, spectrum, insertionOffset);
    }

    return spectrum;
}


std::size_t Spectrum::size() const noexcept
{
    assert(charges.size() == energies.size());
    assert(charges.size() == eigenStateIdxs.size());
    assert(charges.size() == eigenStateCoeffs.size());
    assert(charges.size() == basis.size());
    return charges.size();
}


/*
 * Given eigenstates
 *   |alpha> = sum_x alpha_x |x>
 *   |gamma> = sum_y gamma_y |y>
 * where |x>, |y> are states in spectrum.basis,
 * the matrix elements of operator A are
 *   A^{alpha,gamma} = sum_x sum_y <x| alpha_x A gamma_y |y>
 *                   = sum_x sum_y alpha_x gamma_y <x|A|y>
 *                   = sum_x sum_y alpha_x gamma_y A^{xy}
 */
DSparseMatrix toEigenspaceMatrix(DMatrix const &matrix, Spectrum const &spectrum)
{
    DSparseMatrix res(matrix.rows(), matrix.columns());

    // Iterate over elements of result A^{alpha,gamma}
    for (std::size_t alpha = 0; alpha < spectrum.size(); ++alpha) {
        for (std::size_t gamma = 0; gamma < spectrum.size(); ++gamma) {
            double elem = 0.0;

            // Iterate over basis states
            for (std::size_t x = 0; x < spectrum.eigenStateIdxs[alpha].size(); ++x) {
                for (std::size_t y = 0; y < spectrum.eigenStateIdxs[gamma].size(); ++y) {
                    elem += spectrum.eigenStateCoeffs[alpha][x]
                            * spectrum.eigenStateCoeffs[gamma][y]
                            * matrix(spectrum.eigenStateIdxs[alpha][x],
                                     spectrum.eigenStateIdxs[gamma][y]);
                }
            }

            if (blaze::abs(elem) > 1e-8) {
                res(alpha, gamma) = elem;
            }
        }
    }

    return res;
}
