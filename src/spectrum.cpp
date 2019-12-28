#include "spectrum.hpp"

#include <blaze/math/DynamicVector.h>
#include <blaze/math/lapack/syev.h>

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


    void computeSubSpectrum(SumState const &basis, int const charge,
                            std::vector<std::pair<int, double>> &out)
    {
        SumOperator hamiltonian{ParticleHop{},
                                HoleHop{},
                                SquaredNumberOperator{}};

        auto matrix = toMatrix(hamiltonian, basis);
        blaze::DynamicVector<double> evals(matrix.rows());
        blaze::syev(matrix, evals, 'N', 'U');

        for (double const eval : evals) {
            out.emplace_back(charge, eval);
        }
    }
}


std::vector<std::pair<int, double>> computeSpectrum(SumState basis)
{

    ChargeOperator Q{};

    std::sort(basis.states(), basis.states()+basis.size(),
              [&Q](State const &a, State const &b) {
                  return Q.computeNumber(a) < Q.computeNumber(b);
              });

    EqualChargeIter eci{basis};
    std::vector<std::pair<int, double>> spectrum;
    while (not eci.finished()) {
        auto const [subBasis, charge] = eci.next();
        computeSubSpectrum(subBasis, charge, spectrum);
    }

    return spectrum;
}
