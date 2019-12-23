#ifndef EXACT_HUBBARD_OPERATOR_HPP
#define EXACT_HUBBARD_OPERATOR_HPP


#include <cstdint>
#include <utility>

#include "state.hpp"


constexpr int countPHBefore(State const &state, std::size_t const site) noexcept
{
    assert(site < state.size());
    int count = 0;
    for (std::size_t i = 0; i < site; ++i) {
        count += state.numberOn(i);
    }
    return count;
}


struct ParticleCreator
{
    std::size_t site;

    [[nodiscard]] constexpr std::pair<double, State> apply(State const &state) const noexcept
    {
        if (state.hasParticleOn(site)) {
            return {0, {}};
        }
        State aux{state};
        aux.addParticleOn(site);
        return {countPHBefore(state, site) % 2 == 0 ? +1.0 : -1.0,
                aux};
    }
};


struct ParticleAnnihilator
{
    std::size_t site;

    [[nodiscard]] constexpr std::pair<double, State> apply(State const &state) const noexcept
    {
        if (not state.hasParticleOn(site)) {
            return {0, {}};
        }
        State aux{state};
        aux.removeParticleOn(site);
        return {countPHBefore(state, site) % 2 == 0 ? +1.0 : -1.0,
                aux};
    }
};


#endif //EXACT_HUBBARD_OPERATOR_HPP
