#ifndef EXACT_HUBBARD_STATE_HPP
#define EXACT_HUBBARD_STATE_HPP


#include <array>
#include <cassert>
#include <type_traits>
#include <vector>

#include "config.hpp"


enum class PH : std::uint32_t
{
    n=0b00, p=0b01, h=0b10, ph=0b11
};


constexpr auto underlying(PH const ph) noexcept
{
    return static_cast<std::underlying_type_t<PH>>(ph);
}


constexpr PH operator&(PH const a, PH const b) noexcept
{
    return PH{underlying(a) & underlying(b)};
}


constexpr PH operator|(PH const a, PH const b) noexcept
{
    return PH{underlying(a) | underlying(b)};
}


constexpr PH operator~(PH const ph) noexcept
{
    return PH{~underlying(ph)};
}


class State
{
    std::array<PH, NSITES> sites_{};


public:
    constexpr bool operator==(State const &other) const noexcept
    {
        for (std::size_t i = 0; i < size(); ++i) {
            if (sites_[i] != other.sites_[i]) {
                return false;
            }
        }
        return true;
    }


    [[nodiscard]] constexpr std::size_t size() const noexcept
    {
        return sites_.size();
    }


    constexpr PH &operator[](std::size_t const site) noexcept
    {
        assert(site < sites_.size());
        return sites_[site];
    }


    constexpr PH operator[](std::size_t const site) const noexcept
    {
        assert(site < sites_.size());
        return sites_[site];
    }


    [[nodiscard]] constexpr bool hasParticleOn(std::size_t const site) const noexcept
    {
        assert(site < sites_.size());
        return underlying(sites_[site] & PH::p) != 0;
    }


    [[nodiscard]] constexpr bool hasHoleOn(std::size_t const site) const noexcept
    {
        assert(site < sites_.size());
        return underlying(sites_[site] & PH::h) != 0;
    }


    [[nodiscard]] constexpr int numberOn(std::size_t const site) const noexcept
    {
        assert(site < sites_.size());
#if defined __GNUG__
        static_assert(std::is_same_v<std::underlying_type_t<PH>, unsigned int>);
        return __builtin_popcount(underlying(sites_[site]));
#else
        return static_cast<int>(hasParticleOn(site)) + static_cast<int>(hasHoleOn(site));
#endif
    }


    constexpr void addParticleOn(std::size_t const site) noexcept
    {
        assert(site < sites_.size());
        sites_[site] = sites_[site] | PH::p;
    }


    constexpr void removeParticleOn(std::size_t const site) noexcept
    {
        assert(site < sites_.size());
        sites_[site] = sites_[site] &~ PH::p;
    }


    constexpr void addHoleOn(std::size_t const site) noexcept {
        assert(site < sites_.size());
        sites_[site] = sites_[site] | PH::h;
    }


    constexpr void removeHoleOn(std::size_t const site) noexcept {
        assert(site < sites_.size());
        sites_[site] = sites_[site] &~ PH::h;
    }
};


constexpr std::size_t size(State const &state) noexcept
{
    return state.size();
}


constexpr double dot(State const &a, State const &b) noexcept
{
    // manual because std::mismatch is not constexpr before c++20
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            return 0.0;
        }
    }
    return 1.0;
}


std::vector<State> fockspaceBasis();

#endif //EXACT_HUBBARD_STATE_HPP
