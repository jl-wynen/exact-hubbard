#ifndef EXACT_HUBBARD_SPECTRUM_HPP
#define EXACT_HUBBARD_SPECTRUM_HPP


#include <utility>
#include <vector>

#include "state.hpp"

std::vector<std::pair<int, double>> computeSpectrum(SumState basis);

#endif //EXACT_HUBBARD_SPECTRUM_HPP
