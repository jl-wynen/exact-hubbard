#ifndef EXACT_HUBBARD_SPECTRUM_HPP
#define EXACT_HUBBARD_SPECTRUM_HPP


#include <utility>
#include <vector>

#include "linalg.hpp"
#include "state.hpp"


struct Spectrum
{
    IVector charges;
    DVector energies;
    DSparseMatrix eigenVectors;


    explicit Spectrum(std::size_t size);


    static Spectrum compute(SumState basis);


    std::size_t size() const noexcept;
};

#endif //EXACT_HUBBARD_SPECTRUM_HPP
