#include "state.hpp"

#include <optional>


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
