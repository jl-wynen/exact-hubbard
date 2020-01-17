#ifndef EXACT_HUBBARD_IO_HPP
#define EXACT_HUBBARD_IO_HPP


#include <filesystem>
#include <ostream>

#include "correlators.hpp"
#include "spectrum.hpp"
#include "state.hpp"

namespace fs = std::filesystem;


/// Output a single PH value.
std::ostream &operator<<(std::ostream &os, PH ph);


/// Output a state.
std::ostream &operator<<(std::ostream &os, State const &state);


/// Output a SumState formatted as a sum with coefficients.
std::ostream &operator<<(std::ostream &os, SumState const &states);


/// Write a Spectrum to file.
void saveSpectrum(fs::path const &fname, Spectrum const &spectrum);


/// Write correlators to file.
void saveCorrelators(fs::path const &fname, Correlators const &correlators);

#endif //EXACT_HUBBARD_IO_HPP
