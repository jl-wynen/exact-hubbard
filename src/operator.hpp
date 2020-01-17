#ifndef EXACT_HUBBARD_OPERATOR_HPP
#define EXACT_HUBBARD_OPERATOR_HPP

/** \file
 * \brief Define physical operators.
 *
 * Every operator should derive from the base template `Operator` using CRTP.
 * The derived class must implement
 *   `void apply_implSingleOutparam(State const &state, SumState &out)`
 * which takes a single state, applies the operator, and appends the output to `out`.
 */

#include <cstdint>
#include <utility>

#include "always_false.hpp"
#include "linalg.hpp"
#include "state.hpp"


/// Count the number of particles *and* holes before and excluding the given site.
constexpr int countPHBefore(State const &state, std::size_t const site) noexcept
{
    assert(site < state.size());
    int count = 0;
    for (std::size_t i = 0; i < site; ++i) {
        count += state.numberOn(i);
    }
    return count;
}


/// Check whether a type has a member function `apply_implSingleOutparam`.
template <typename T, typename = void>
struct hasApplyImplSingleOutparam : std::false_type {};

template <typename T>
struct hasApplyImplSingleOutparam<T,
        std::void_t<decltype(std::declval<T>().apply_implSingleOutparam(std::declval<State>(),
                                                                        std::declval<SumState&>()))>>
        : std::true_type {};


/// Base template for operators.
/**
 * Operators implementations should inherit from this class using CRTP.
 */
template <typename T>
struct Operator
{
    /// The type of the operator implementation.
    using Derived = T;


    /// Convert `*this` to a reference to `Derived`.
    [[nodiscard]] constexpr Derived &asDerived() noexcept
    {
        return static_cast<Derived&>(*this);
    }


    /// Convert `*this` to a reference to `Derived`.
    [[nodiscard]] constexpr Derived const &asDerived() const noexcept
    {
        return static_cast<Derived const&>(*this);
    }


    /// Apply the operator to a single state and append the result to `out`.
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


    /// Apply the operator to all states in `in` and append the results to `out`.
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


    /// Apply the operator to all states in `states` and return a new SumState with the results.
    [[nodiscard]] SumState apply(SumState const &states) const
    {
        SumState out;
        // Might over allocate if fewer states are produced but that should be OK.
        out.reserve(states.size());
        apply(states, out);
        return out;
    }
};


/// The sum of multiple operators.
/**
 * Use e.g. as
 * ```{.cpp}
SumOperator hamiltonian{ParticleHop{},
                        HoleHop{},
                        SquaredNumberOperator{}};
   ```
 * `hamiltonian` can then be applied like any other operator.
 */
template <typename... Operators>
struct SumOperator : Operator<SumOperator<Operators...>>
{
    /// Stores all summands (sub operators).
    std::tuple<Operator<Operators>...> operators;


    /// Construct from one or more summands.
    explicit constexpr SumOperator(Operator<Operators> const & ... ops)
            : operators{ops...}
    { }


    /// Implementation of apply.
    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        if constexpr (sizeof...(Operators) > 0) {
            doApply<sizeof...(Operators)-1>(state, out);
        }
    }


private:
    /// Recursively apply sub operators.
    template <std::size_t Idx>
    void doApply(State const &state, SumState &out) const
    {
        std::get<Idx>(operators).apply(state, out);

        if constexpr (Idx > 0) {
            doApply<Idx-1>(state, out);
        }
    }
};


/// Deduction guide for SumOperator.
template <typename... Operators>
SumOperator(Operator<Operators> const & ...) -> SumOperator<Operators...>;


/// Creator for a single particle at a given site.
struct ParticleCreator : Operator<ParticleCreator>
{
    /// Lattice site the creator operates on.
    std::size_t site;


    /// Specify the site at which particles are created.
    explicit constexpr ParticleCreator(std::size_t const s) noexcept : site{s} { }


    /// Implementation of apply.
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


/// Annihilator for a single particle at a given site.
struct ParticleAnnihilator : Operator<ParticleAnnihilator>
{
    /// Lattice site the annihilator operates on.
    std::size_t site;


    /// Specify the site at which particles are annihilated.
    explicit constexpr ParticleAnnihilator(std::size_t const s) noexcept : site{s} { }


    /// Implementation of apply.
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


/// Creator for a single hole at a given site.
struct HoleCreator : Operator<HoleCreator>
{
    /// Lattice site the creator operates on.
    std::size_t site;


    /// Specify the site at which particles are created.
    explicit constexpr HoleCreator(std::size_t const s) noexcept : site{s} { }


    /// Implementation of apply.
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


/// Annihilator for a single hole at a given site.
struct HoleAnnihilator : Operator<HoleAnnihilator>
{
    /// Lattice site the annihilator operates on.
    std::size_t site;


    /// Specify the site at which holes are annihilated.
    explicit constexpr HoleAnnihilator(std::size_t const s) noexcept : site{s} { }


    /// Implementation of apply.
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
 * \f$ U/2 \sum_x (n_x - \tilde{n}_x)^2 \f$
 * Where \f$ n_x \f$ and \f$ \tilde{n} \f$ count particles and holes at site x.
 * \tparam USE_PREFACTOR If `true`, multiply result by U/2 as shown in equation above.
 */
template <bool USE_PREFACTOR=true>
struct SquaredNumberOperator : Operator<SquaredNumberOperator<USE_PREFACTOR>>
{
    /// Implementation of apply.
    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        int number = 0;
        for (std::size_t site = 0; site < NSITES; ++site) {
            if ((state.hasParticleOn(site) ^ state.hasHoleOn(site)) != 0) {
                number++;
            }
        }
        if (number != 0) {
            if constexpr (USE_PREFACTOR) {
                out.push(U / 2.0 * static_cast<double>(number), state);
            }
            else {
                out.push(static_cast<double>(number), state);
            }
        }
    }
};


/**
 * \f$ \sum_x (n_x - \tilde{n}_x) \f$
 */
struct ChargeOperator : Operator<ChargeOperator>
{
    /// Implementation of apply.
    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        int const number = computeCharge(state);
        if (number != 0) {
            out.push(static_cast<double>(number), state);
        }
    }


    /// Compute the charge of a state.
    [[nodiscard]] int computeCharge(State const &state) const noexcept
    {
        int number = 0;
        for (std::size_t site = 0; site < NSITES; ++site) {
            number += static_cast<int>(state.hasParticleOn(site))
                      - static_cast<int>(state.hasHoleOn(site));
        }
        return number;
    }
};


/**
 * Hopping operator for particles.
 * Annihilates and creates particles according to the hopping matrix
 * as parameterised by nearestNeighbours and kappa.
 * The operator can be written as
 * \f[
   -\kappa \sum_{\langle x, y\rangle}\, a_x^\dagger a_y
 * \f]
 */
struct ParticleHop : Operator<ParticleHop>
{
    using Operator<ParticleHop>::apply;


    /// Implementation of apply.
    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        for (auto const [a, b] : nearestNeighbours) {
            if (state.hasParticleOn(a) and not state.hasParticleOn(b)) {
                auto const [coef, newState] = doHop(state, a, b);
                out.push(coef, newState);
            }
            // Can use else here because we can never hop to *and* from a site.
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


/**
 * Hopping operator for particles.
 * Annihilates and creates particles according to the hopping matrix
 * as parameterised by nearestNeighbours and kappa.
 * The operator can be written as
 * \f[
   \kappa \sum_{\langle x, y\rangle}\, b_x^\dagger b_y
 * \f]
 */
struct HoleHop : Operator<HoleHop>
{
    using Operator<HoleHop>::apply;


    /// Implementation of apply.
    void apply_implSingleOutparam(State const &state, SumState &out) const
    {
        for (auto const [a, b] : nearestNeighbours) {
            if (state.hasHoleOn(a) and not state.hasHoleOn(b)) {
                auto const [coef, newState] = doHop(state, a, b);
                out.push(coef, newState);
            }
                // Can use else here because we can never hop to *and* from a site.
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

        return {kappa * ((nSwapAnnihilate+nSwapCreate) % 2 == 0 ? +1.0 : -1.0),
                newState};
    }
};


/**
 * Compute all matrix elements of an operator.

 * @param op %Operator \f$ O \f$.
 * @param basis Each state in `basis` is a basis state \f$ |i\rangle \f$
 * @return \f$ M_{ij} = \langle i | O | j \rangle \f$
 */
template <typename T>
DMatrix toMatrix(Operator<T> const &op, SumState const &basis)
{
    DMatrix mat{basis.size(), basis.size()};
    SumState out;

    for (std::size_t j = 0; j < mat.columns(); ++j)
    {
        out.clear();
        // |j>
        auto const &[coefj, statej] = basis[j];
        // out = op |j>
        op.apply(statej, out);

        for (std::size_t i = 0; i < mat.rows(); ++i) {
            // |i>
            auto const &[coefi, statei] = basis[i];
            // <i|op|j>
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
