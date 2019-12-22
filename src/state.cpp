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


std::vector<State> fockspaceBasis()
{
    State auxState;
    std::vector<State> states{auxState};

    while (true) {
        if (not increment(auxState)) {
            break;
        }
        states.emplace_back(auxState);
    }

    return states;
}
