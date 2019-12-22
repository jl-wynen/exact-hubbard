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


std::ostream &operator<<(std::ostream &os, State const &state) {
    os << '|';
    for (std::size_t i = 0; i + 1 < state.size(); ++i) {
        os << state[i] << ' ';
    }
    os << state[state.size()-1] << '>';
    return os;
}

#endif //EXACT_HUBBARD_IO_HPP
