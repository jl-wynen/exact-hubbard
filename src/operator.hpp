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


template <typename T, typename = void>
struct hasSingleStateApplyWithOut : std::false_type {};

template <typename T>
struct hasSingleStateApplyWithOut<T,
        std::void_t<decltype(std::declval<T>().apply(std::declval<SumState>(), std::declval<SumState&>()))>>
        : std::true_type {};


template <typename T>
struct Operator
{
    using Derived = T;


    [[nodiscard]] constexpr Derived &asDerived() noexcept
    {
        return static_cast<Derived&>(*this);
    }


    [[nodiscard]] constexpr Derived const &asDerived() const noexcept
    {
        return static_cast<Derived const&>(*this);
    }


    void apply(State const &state, SumState &out) const
    {
        if constexpr (hasSingleStateApplyWithOut<Derived>::value) {
            // use implementation provided by derived type
            asDerived().apply(state, out);
        }
        else {
            // use generic implementation
            auto const[c, s] = asDerived().apply(state);
            if (c != 0.0) {
                out.push(c, s);
            }
        }
    }


    void apply(SumState const &in, SumState &out) const
    {
        std::size_t outa = out.size();
        for (std::size_t i = 0; i < in.size(); ++i) {
            auto const &[coef, state] = in[i];
            apply(state, out);

            for (std::size_t const outb = out.size(); outa < outb; ++outa) {
                out[outa].first *= coef;
            }
        }
    }


    [[nodiscard]] SumState apply(SumState const &states) const
    {
        SumState out;
        // Might over allocate if fewer states are produced but that should be OK.
        out.reserve(states.size());
        apply(states, out);
        return out;
    }
};


struct ParticleCreator : Operator<ParticleCreator>
{
    std::size_t site;


    explicit constexpr ParticleCreator(std::size_t const s) noexcept : site{s} { }


    using Operator<ParticleCreator>::apply;


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


struct ParticleAnnihilator : Operator<ParticleAnnihilator>
{
    std::size_t site;


    explicit constexpr ParticleAnnihilator(std::size_t const s) noexcept : site{s} { }


    using Operator<ParticleAnnihilator>::apply;


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


/**
 * (n_x - tilde{n}_x)^2
 */
struct NumberOperator : Operator<NumberOperator>
{
    std::size_t site;

    explicit constexpr NumberOperator(std::size_t const s) noexcept : site{s} { }


    using Operator<NumberOperator>::apply;


    [[nodiscard]] constexpr std::pair<double, State> apply(State const &state) const noexcept
    {
        return {(state.hasParticleOn(site) ^ state.hasHoleOn(site)) ? 1.0 : 0.0,
                state};
    }
};


struct ParticleHop : Operator<ParticleHop>
{
    using Operator<ParticleHop>::apply;


    void apply(State const &state, SumState &out) const
    {
        for (auto const [a, b] : nearestNeighbours) {
            if (state.hasParticleOn(a) and not state.hasParticleOn(b)) {
                auto const [coef, newState] = doHop(state, a, b);
                out.push(coef, newState);
            }
            // Can use else here because we can never hop to _and_ from a site.
            else if (state.hasParticleOn(b) and not state.hasParticleOn(a)) {
                auto const [coef, newState] = doHop(state, b, a);
                out.push(coef, newState);
            }
        }
    }


private:
    [[nodiscard]] std::pair<double, State>
            doHop(State const &state, std::size_t const from, std::size_t const to) const noexcept
    {
        // How often does the annihilator have to swap places with another operator?
        auto const nSwapAnnihilate = countPHBefore(state, from);
        // Annihilate
        State newState{state};
        newState.removeParticleOn(from);

        // How often does the creator have to swap places with another operator?
        auto const nSwapCreate = countPHBefore(newState, to);
        // Create
        newState.addParticleOn(to);

        return {-kappa * ((nSwapAnnihilate+nSwapCreate) % 2 == 0 ? +1.0 : -1.0),
                newState};
    }
};

#endif //EXACT_HUBBARD_OPERATOR_HPP
