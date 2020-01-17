#ifndef EXACT_HUBBARD_STATE_HPP
#define EXACT_HUBBARD_STATE_HPP

/** \file
 * \brief Classes State and SumState.
 */

#include <array>
#include <cassert>
#include <type_traits>
#include <vector>

#include "config.hpp"


/// Indicate the presence of particles and holes at a lattice site.
enum class PH : std::uint32_t
{
    n=0b00,  ///< Neither particle or hole
    p=0b01,  ///< A particle, no hole
    h=0b10,  ///< A hole, no particle
    ph=0b11  ///< Both particle and hole
};


/// Case a PH to its underlying type.
constexpr auto underlying(PH const ph) noexcept
{
    return static_cast<std::underlying_type_t<PH>>(ph);
}


/// Bit-wise and on PH.
constexpr PH operator&(PH const a, PH const b) noexcept
{
    return PH{underlying(a) & underlying(b)};
}


/// Bit-wise or on PH.
constexpr PH operator|(PH const a, PH const b) noexcept
{
    return PH{underlying(a) | underlying(b)};
}


/// Bit-wise 2's complement on PH.
constexpr PH operator~(PH const ph) noexcept
{
    return PH{~underlying(ph)};
}


/**
 * Store PH values for all NSITES lattice sites.
 */
class State
{
    std::array<PH, NSITES> sites_{};


public:
    /// Return `true` iff *all* sites are equal.
    constexpr bool operator==(State const &other) const noexcept
    {
        for (std::size_t i = 0; i < size(); ++i) {
            if (sites_[i] != other.sites_[i]) {
                return false;
            }
        }
        return true;
    }


    /// Return number of alttice sites.
    [[nodiscard]] constexpr std::size_t size() const noexcept
    {
        return sites_.size();
    }


    /// Access PH value at a site.
    constexpr PH &operator[](std::size_t const site) noexcept
    {
        assert(site < sites_.size());
        return sites_[site];
    }


    /// Access PH value at a site.
    constexpr PH operator[](std::size_t const site) const noexcept
    {
        assert(site < sites_.size());
        return sites_[site];
    }


    /// Return `true` if there is a particle on a site.
    [[nodiscard]] constexpr bool hasParticleOn(std::size_t const site) const noexcept
    {
        assert(site < sites_.size());
        return underlying(sites_[site] & PH::p) != 0;
    }


    /// Return `true` if there is a heol on a site.
    [[nodiscard]] constexpr bool hasHoleOn(std::size_t const site) const noexcept
    {
        assert(site < sites_.size());
        return underlying(sites_[site] & PH::h) != 0;
    }


    /// Return the number of particles+holes on a site.
    [[nodiscard]] constexpr int numberOn(std::size_t const site) const noexcept
    {
        assert(site < sites_.size());
#if defined __GNUG__
        // count number of set bits if possible
        static_assert(std::is_same_v<std::underlying_type_t<PH>, unsigned int>);
        return __builtin_popcount(underlying(sites_[site]));
#else
        // count manually (does not get optimised to popcount automatically)
        return static_cast<int>(hasParticleOn(site)) + static_cast<int>(hasHoleOn(site));
#endif
    }


    /// Make it so there is a particle at a site regardless of whether there was one already.
    constexpr void addParticleOn(std::size_t const site) noexcept
    {
        assert(site < sites_.size());
        sites_[site] = sites_[site] | PH::p;
    }


    /// Make it so there is no particle at a site regardless of whether there was one in the first place.
    constexpr void removeParticleOn(std::size_t const site) noexcept
    {
        assert(site < sites_.size());
        sites_[site] = sites_[site] &~ PH::p;
    }


    /// Make it so there is a hole at a site regardless of whether there was one already.
    constexpr void addHoleOn(std::size_t const site) noexcept {
        assert(site < sites_.size());
        sites_[site] = sites_[site] | PH::h;
    }


    /// Make it so there is no hole at a site regardless of whether there was one in the first place.
    constexpr void removeHoleOn(std::size_t const site) noexcept {
        assert(site < sites_.size());
        sites_[site] = sites_[site] &~ PH::h;
    }
};


/// Return the number of sites based on a state.
constexpr std::size_t size(State const &state) noexcept
{
    return state.size();
}


/// Compute the dot product of two states.
constexpr double dot(State const &a, State const &b) noexcept
{
    return a == b ? +1.0 : 0.0;
}


/**
 * Store multiple states and coefficients and represent them as their sum.
 *
 * Corresponds to sum_i sumState[i][0]*sumState[i][1].
 */
class SumState
{
    std::vector<double> coefs_;
    std::vector<State> states_;

public:
    /// Reserve memory for n states.
    void reserve(std::size_t const n)
    {
        coefs_.reserve(n);
        states_.reserve(n);
    }


    /// Return the number of currently stored states.
    [[nodiscard]] std::size_t size() const noexcept
    {
        return states_.size();
    }


    /// Return coefficient and state number `i`.
    std::pair<double&, State&> operator[](std::size_t const i) noexcept
    {
        return {coefs_[i], states_[i]};
    }


    /// Return coefficient and state number `i`.
    std::pair<double, State const&> operator[](std::size_t const i) const noexcept
    {
        return {coefs_[i], states_[i]};
    }


    /// Access underlying state storage.
    [[nodiscard]] State const *states() const noexcept
    {
        return states_.data();
    }


    /// Access underlying state storage.
    [[nodiscard]] State *states() noexcept
    {
        return states_.data();
    }


    /// Append a new state.
    void push(double const coef, State const &state)
    {
        states_.emplace_back(state);
        coefs_.push_back(coef);
    }


    /// Remove all stored states and coefficients.
    void clear() noexcept
    {
        coefs_.clear();
        states_.clear();
    }


    /// Remove duplicate states by adding up their coefficients.
    void compress();
};


/// Compute the dot product of two SumStates.
inline double dot(SumState const &a, SumState const &b) noexcept
{
    double res = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        auto const &[ca, sa] = a[i];
        for (std::size_t j = 0; j < b.size(); ++j) {
            auto const &[cb, sb] = b[j];
            res += ca * cb * dot(sa, sb);
        }
    }
    return res;
}


/**
 * Construct the basis of the fockspace.
 * \return SumState containing basis vectors.
 *         All coefficients are `1.0`.
 */
SumState fockspaceBasis();

#endif //EXACT_HUBBARD_STATE_HPP
