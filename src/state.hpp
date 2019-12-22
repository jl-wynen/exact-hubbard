#ifndef EXACT_HUBBARD_STATE_HPP
#define EXACT_HUBBARD_STATE_HPP


#include <array>
#include <cassert>

#include "config.hpp"


enum class PH : int
{
    n=0b00, p=0b01, h=0b10, ph=0b11
};


struct State
{
    std::array<PH, NSITES> sites{};


    constexpr PH &operator[](size_t const i) noexcept
    {
        assert(i < sites.size());
        return sites[i];
    }


    constexpr PH operator[](size_t const i) const noexcept
    {
        assert(i < sites.size());
        return sites[i];
    }
};

#endif //EXACT_HUBBARD_STATE_HPP
