#include <iostream>
#include <fstream>
#include <optional>

#include "state.hpp"
#include "spectrum.hpp"


void computeNonInteractingSpectrum()
{
    static_assert(U == 0.0);

    auto const fockspace = fockspaceBasis();
    auto const spectrum = computeSpectrum(fockspace);
    std::ofstream ofs("../spectrum.dat");
    ofs << "#  Q  E\n";
    for (auto const [charge, energy] : spectrum) {
        ofs << charge << ' ' << energy << '\n';
    }
}

int main()
{
    computeNonInteractingSpectrum();
}
