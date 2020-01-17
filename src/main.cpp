#include <iostream>
#include <chrono>

#include "correlators.hpp"
#include "io.hpp"
#include "spectrum.hpp"


int main()
{
    std::cout << "Nx = " << NSITES << ",  Nt = " << NT << '\n'
              << "beta = " << beta << ",  U = " << U << ",  kappa = " << kappa << '\n';

    // spectrum
    auto startTime = std::chrono::high_resolution_clock::now();
    auto const spectrum = Spectrum::compute(fockspaceBasis());
    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "Time to compute spectrum: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                      endTime-startTime
              ).count() << "ms\n";
    saveSpectrum("../spectrum.dat", spectrum);

    // correlators
    startTime = std::chrono::high_resolution_clock::now();
    auto const correlators = computeCorrelators(spectrum);
    endTime = std::chrono::high_resolution_clock::now();
    std::cout << "Time to compute correlators: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                      endTime-startTime
              ).count() << "ms\n";
    saveCorrelators("../correlators.dat", correlators);
}
