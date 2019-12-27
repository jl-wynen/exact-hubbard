#ifndef EXACT_HUBBARD_OPERATOR_HPP
#define EXACT_HUBBARD_OPERATOR_HPP


#include <cstdint>
#include <utility>

#include "always_false.hpp"
#include "matrix.hpp"
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
struct hasApplyImplSingleOutparam : std::false_type {};

template <typename T>
struct hasApplyImplSingleOutparam<T,
        std::void_t<decltype(std::declval<T>().apply_implSingleOutparam(std::declval<State>(),
                                                                        std::declval<SumState&>()))>>
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
        if constexpr (hasApplyImplSingleOutparam<Derived>::value) {
            asDerived().apply_implSingleOutparam(state, out);
        }
        else {
            static_assert(alwaysFalse_v<Derived>,
                          "The operator class has no valid implementation of apply");
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


template <typename... Operators>
struct SumOperator : Operator<SumOperator<Operators...>>
{
    std::tuple<Operator<Operators>...> operators;


    explicit constexpr SumOperator(Operator<Operators> const & ... ops)
            : operators{ops...}
    { }


    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        if constexpr (sizeof...(Operators) > 0) {
            doApply<sizeof...(Operators)-1>(state, out);
        }
    }


private:
    template <std::size_t Idx>
    void doApply(State const &state, SumState &out) const
    {
        std::get<Idx>(operators).apply(state, out);

        if constexpr (Idx > 0) {
            doApply<Idx-1>(state, out);
        }
    }
};

template <typename... Operators>
SumOperator(Operator<Operators> const & ...) -> SumOperator<Operators...>;


struct ParticleCreator : Operator<ParticleCreator>
{
    std::size_t site;


    explicit constexpr ParticleCreator(std::size_t const s) noexcept : site{s} { }


    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        if (state.hasParticleOn(site)) {
            return;  // cannot create where there is a particle already
        }
        State aux{state};
        aux.addParticleOn(site);
        out.push(countPHBefore(state, site) % 2 == 0 ? +1.0 : -1.0,
                 aux);
    }
};


struct ParticleAnnihilator : Operator<ParticleAnnihilator>
{
    std::size_t site;


    explicit constexpr ParticleAnnihilator(std::size_t const s) noexcept : site{s} { }


    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        if (not state.hasParticleOn(site)) {
            return;  // cannot destroy a particle when there is none
        }
        State aux{state};
        aux.removeParticleOn(site);
        out.push(countPHBefore(state, site) % 2 == 0 ? +1.0 : -1.0,
                 aux);
    }
};


struct HoleCreator : Operator<HoleCreator>
{
    std::size_t site;


    explicit constexpr HoleCreator(std::size_t const s) noexcept : site{s} { }


    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        if (state.hasHoleOn(site)) {
            return;  // cannot create a hole where there is a hole already
        }
        State aux{state};
        aux.addHoleOn(site);
        auto const nswaps = countPHBefore(state, site) + (state.hasParticleOn(site) ? 1 : 0);
        out.push(nswaps % 2 == 0 ? +1.0 : -1.0,
                 aux);
    }
};


struct HoleAnnihilator : Operator<HoleAnnihilator>
{
    std::size_t site;


    explicit constexpr HoleAnnihilator(std::size_t const s) noexcept : site{s} { }


    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        if (not state.hasHoleOn(site)) {
            return;  // cannot destroy a hole when there is none
        }
        State aux{state};
        aux.removeHoleOn(site);
        auto const nswaps = countPHBefore(state, site) + (state.hasParticleOn(site) ? 1 : 0);
        out.push(nswaps % 2 == 0 ? +1.0 : -1.0,
                 aux);
    }
};


/**
 * U/2 sum_x (n_x - tilde{n}_x)^2
 */
struct GlobalNumberOperator : Operator<GlobalNumberOperator>
{
    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        int number = 0;
        for (std::size_t site = 0; site < NSITES; ++site) {
            if ((state.hasParticleOn(site) ^ state.hasHoleOn(site)) != 0) {
                number++;
            }
        }
        if (number != 0) {
            out.push(U / 2.0 * static_cast<double>(number), state);
        }
    }
};


struct ParticleHop : Operator<ParticleHop>
{
    using Operator<ParticleHop>::apply;


    void apply_implSingleOutparam(State const &state, SumState &out) const
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


struct HoleHop : Operator<HoleHop>
{
    using Operator<HoleHop>::apply;


    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        for (auto const [a, b] : nearestNeighbours) {
            if (state.hasHoleOn(a) and not state.hasHoleOn(b)) {
                auto const [coef, newState] = doHop(state, a, b);
                out.push(coef, newState);
            }
                // Can use else here because we can never hop to _and_ from a site.
            else if (state.hasHoleOn(b) and not state.hasHoleOn(a)) {
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
        auto const nSwapAnnihilate = countPHBefore(state, from) + (state.hasParticleOn(from) ? 1 : 0);
        // Annihilate
        State newState{state};
        newState.removeHoleOn(from);

        // How often does the creator have to swap places with another operator?
        auto const nSwapCreate = countPHBefore(newState, to) + (state.hasParticleOn(to) ? 1 : 0);
        // Create
        newState.addHoleOn(to);

        return {-kappa * ((nSwapAnnihilate+nSwapCreate) % 2 == 0 ? +1.0 : -1.0),
                newState};
    }
};


/**
 * M_{ij} = <i|O|j>
 * @tparam T
 * @param op
 * @param basis
 * @return
 */
template <typename T>
Matrix<double> toMatrix(Operator<T> const &op, SumState const &basis)
{
    Matrix<double> mat{basis.size(), basis.size()};
    SumState out;

    for (std::size_t j = 0; j < mat.ncol(); ++j)
    {
        out.clear();
        auto const &[coefj, statej] = basis[j];
        op.apply(statej, out);

        for (std::size_t i = 0; i < mat.nrow(); ++i) {
            auto const &[coefi, statei] = basis[i];
            double matelem = 0.0;
            for (std::size_t k = 0; k < out.size(); ++k) {
                auto const &[coefk, statek] = out[k];
                matelem += dot(statek, statei) * coefk * coefi * coefj;
            }
            mat(i, j) = matelem;
        }
    }

    return mat;
}

#endif //EXACT_HUBBARD_OPERATOR_HPP
