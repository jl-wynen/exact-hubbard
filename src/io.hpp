#ifndef EXACT_HUBBARD_IO_HPP
#define EXACT_HUBBARD_IO_HPP


#include <ostream>

#include "state.hpp"


std::ostream &operator<<(std::ostream &os, PH ph);


std::ostream &operator<<(std::ostream &os, State const &state);


std::ostream &operator<<(std::ostream &os, SumState const &states);

#endif //EXACT_HUBBARD_IO_HPP
