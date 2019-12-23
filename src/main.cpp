#include <iostream>
#include <iomanip>

#include "state.hpp"
#include "io.hpp"
#include "operator.hpp"


int main()
{
    auto const fockspace = fockspaceBasis();

    ParticleCreator ad{1};
    ParticleAnnihilator a{1};


    for (size_t i = 0; i < fockspace.size(); ++i) {
        auto [c, s] = fockspace[i];
        std::cout << std::setw(2) << c << ' ' <<  s << '\n';
    }
    std::cout << fockspace.size() << '\n';

}
