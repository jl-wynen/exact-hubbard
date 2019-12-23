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
        auto const [c, s] = asDerived().apply(state);
        if (c != 0.0) {
            out.push(c, s);
        }
    }


    void apply(SumState const &in, SumState &out) const
    {
        std::size_t outa = out.size();
        for (std::size_t i = 0; i < in.size(); ++i) {
            auto const &[coef, state] = in[i];
            apply(state, out);

            std::size_t outb = out.size();
            for (; outa < outb; ++outa) {
                out[outa].first *= coef;
            }
        }
    }


    [[nodiscard]] SumState apply(SumState const &states) const
    {
        SumState out;
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


struct NumberOperator : Operator<NumberOperator>
{
    std::size_t site;

    explicit constexpr NumberOperator(std::size_t const s) noexcept : site{s} { }


    using Operator<NumberOperator>::apply;


    [[nodiscard]] constexpr std::pair<double, State> apply(State const &state) const noexcept
    {
        return {static_cast<double>(state.numberOn(site)), state};
    }
};

#endif //EXACT_HUBBARD_OPERATOR_HPP
