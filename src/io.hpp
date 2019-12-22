#ifndef EXACT_HUBBARD_IO_HPP
#define EXACT_HUBBARD_IO_HPP


#include <ostream>

#include "state.hpp"


std::ostream &operator<<(std::ostream &os, PH ph)
{
    switch (ph) {
        case PH::n:
            os << ". "; break;
        case PH::p:
            os << "p "; break;
        case PH::h:
            os << " h"; break;
        case PH::ph:
            os << "ph"; break;
    }
    return os;
}


std::ostream &operator<<(std::ostream &os, State const &state)
{
    os << '|';
    for (auto x : state.sites) {
        os << x << ' ';
    }
    os << '>';
    return os;
}

#endif //EXACT_HUBBARD_IO_HPP
