#include <iostream>
#include <iomanip>

#include "state.hpp"
#include "io.hpp"
#include "operator.hpp"


int main()
{
    auto const fockspace = fockspaceBasis();

    ParticleCreator ad{1};
    ParticleAnnihilator a{0};

    auto aux = ad.apply(fockspace);
    auto out = a.apply(aux);

    for (size_t i = 0; i < fockspace.size(); ++i) {
        auto [c, s] = fockspace[i];
        auto [ca, sa] = aux[i];
        auto [c2, s2] = out[i];
        std::cout << std::setw(2) << c << ' ' <<  s
                  << "   " << std::setw(2) << ca << ' ' << sa
                  << "   " << std::setw(2) << c2 << ' ' << s2 << '\n';
    }
    std::cout << fockspace.size() << '\n';

}
