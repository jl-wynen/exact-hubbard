#include "state.hpp"

#include <optional>


namespace {
    template <typename T>
    void erase(std::vector<T> &vec, std::size_t const i)
    {
        vec.erase(vec.cbegin()
                  + static_cast<typename decltype(vec.cbegin())::difference_type>(i));
    }
}


void SumState::compress()
{
    for (std::size_t i = 0; i < coefs_.size(); ++i) {
        // join elements with the same state
        for (std::size_t j = coefs_.size()-1; j > i; --j) {
            if (states_[i] == states_[j]) {
                coefs_[i] += coefs_[j];
                erase(coefs_, j);
                erase(states_, j);
            }
        }

        // erase elements with coefficient 0
        if (std::abs(coefs_[i]) < 1e-13) {
            erase(coefs_, i);
            erase(states_, i);
            --i;  // check the new element at this position next
        }
    }
}

namespace {
    constexpr std::optional<PH> increment(PH ph) noexcept
    {
        if (ph == PH::ph) {
            return {};
        }

        return PH{underlying(ph) + 1};
    }


    [[nodiscard]] constexpr bool increment(State &state, size_t const dim = 0) noexcept
    {
        auto const inc = increment(state[dim]);
        if (inc.has_value()) {
            state[dim] = *inc;
            return true;
        }
        else {
            state[dim] = PH::n;
            if (dim == state.size()-1) {
                return false;
            }
            return increment(state, dim+1);
        }
    }
}


SumState fockspaceBasis()
{
    SumState basis;
    State auxState;
    basis.push(1.0, auxState);

    while (true) {
        if (not increment(auxState)) {
            break;
        }
        basis.push(1.0, auxState);
    }

    return basis;
}
