#include "io.hpp"

#include <fstream>



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


std::ostream &operator<<(std::ostream &os, State const &state)
{
    os << '|';
    for (std::size_t i = 0; i + 1 < state.size(); ++i) {
        os << state[i] << ' ';
    }
    os << state[state.size()-1] << '>';
    return os;
}


std::ostream &operator<<(std::ostream &os, SumState const &states)
{
    for (std::size_t i = 0; i < states.size(); ++i) {
        auto const [coef, state] = states[i];
        os << coef << state;
        if (i+1 < states.size()) {
            os << " + ";
        }
    }
    return os;
}


void saveSpectrum(fs::path const &fname, Spectrum const &spectrum)
{
    std::ofstream ofs(fname);
    ofs << "#  Q  E\n";
    for (std::size_t i = 0; i < spectrum.size(); ++i) {
        ofs << spectrum.charges[i] << ' ' << spectrum.energies[i] << '\n';
    }
}


void saveCorrelators(fs::path const &fname, Correlators const &correlators)
{
    std::ofstream ofs{fname};
    ofs << "#~ correlator\n#  nx  nt\n"
        << NSITES << ' ' << NT
        << "\n#  U  kappa  beta\n"
        << U << ' ' << kappa << ' ' << beta
        << "\n#  data\n";
    for (auto const x : correlators.data) {
        ofs << x << ' ';
    }
}
