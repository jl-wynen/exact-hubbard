#ifndef EXACT_HUBBARD_STATE_HPP
#define EXACT_HUBBARD_STATE_HPP


#include <array>
#include <cassert>

#include "config.hpp"


enum class PH : int
{
    n=0b00, p=0b01, h=0b10, ph=0b11
};


class State
{
    std::array<PH, NSITES> sites_{};


public:
    [[nodiscard]] constexpr std::size_t size() const noexcept
    {
        return sites_.size();
    }


    constexpr PH &operator[](size_t const i) noexcept
    {
        assert(i < sites_.size());
        return sites_[i];
    }


    constexpr PH operator[](size_t const i) const noexcept
    {
        assert(i < sites_.size());
        return sites_[i];
    }
};


constexpr std::size_t size(State const &state) noexcept
{
    return state.size();
}


#endif //EXACT_HUBBARD_STATE_HPP
