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
            currentCharge_ = Q_.computeNumber(basis_[currentI_].second);
            res.push(basis_[currentI_].first, basis_[currentI_].second);
            currentI_++;

            for (; currentI_ < basis_.size()
                   and Q_.computeNumber(basis_[currentI_].second) == currentCharge_;
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
            out.energies[insertionOffset + i] = evals[i] / kappa;

            // `syev` stores the eigenvectors row-wise in `matrix`.
//            out.eigenStates[insertionOffset + i] = stateInBasis(blaze::row(matrix, i), basis);
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
                  return Q.computeNumber(a) < Q.computeNumber(b);
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
