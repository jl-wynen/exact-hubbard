#ifndef EXACT_HUBBARD_IO_HPP
#define EXACT_HUBBARD_IO_HPP


#include <ostream>

#include "matrix.hpp"
#include "state.hpp"


std::ostream &operator<<(std::ostream &os, PH ph);


std::ostream &operator<<(std::ostream &os, State const &state);


std::ostream &operator<<(std::ostream &os, DMatrix const &matrix);

#endif //EXACT_HUBBARD_IO_HPP
