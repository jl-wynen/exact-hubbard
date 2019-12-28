#include <iostream>
#include <fstream>
#include <optional>

#include "state.hpp"
#include "io.hpp"
#include "operator.hpp"
#include "linalg.hpp"
#include "spectrum.hpp"


int main()
{
    auto const fockspace = fockspaceBasis();

    auto const spectrum = computeSpectrum(fockspace);
    std::ofstream ofs("../spectrum.dat");
    for (auto const [n, e] : spectrum) {
        ofs << n << ' ' << e << '\n';
    }
}
