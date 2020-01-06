#include <iostream>
#include <fstream>
#include <optional>

#include "state.hpp"
#include "spectrum.hpp"


void computeNonInteractingSpectrum()
{
    static_assert(U == 0.0);

    auto const fockspace = fockspaceBasis();
    auto const spectrum = Spectrum::compute(fockspace);
    std::ofstream ofs("../spectrum.dat");
    ofs << "#  Q  E\n";
    for (std::size_t i = 0; i < spectrum.size(); ++i) {
        ofs << spectrum.charges[i] << ' ' << spectrum.energies[i] << '\n';
    }
}

int main()
{
    computeNonInteractingSpectrum();
}
