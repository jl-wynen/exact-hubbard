#include <iostream>
#include <optional>
#include <vector>

#include "state.hpp"
#include "io.hpp"


constexpr std::optional<PH> increment(PH ph) noexcept
{
    if (ph == PH::ph) {
        return {};
    }

    return PH{static_cast<std::underlying_type_t<PH>>(ph) + 1};
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


std::vector<State> fockspaceBasis()
{
    State auxState;
    std::vector<State> states{auxState};

    while (true) {
        std::cout << auxState << '\n';
        if (not increment(auxState)) {
            break;
        }
        states.emplace_back(auxState);
    }

    return states;
}


int main()
{
    auto fockspace = fockspaceBasis();

    for (auto s : fockspace) {
        std::cout << s << '\n';
    }
    std::cout << fockspace.size() << '\n';

}
