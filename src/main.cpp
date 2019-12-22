#include <iostream>
#include <iomanip>

#include "state.hpp"
#include "io.hpp"
#include "operator.hpp"


int main()
{
    auto fockspace = fockspaceBasis();

    ParticleCreator ad{1};

    for (auto s : fockspace) {
        auto [c, ss] = ad.apply(s);
        std::cout << s << ' ' << std::setw(2) << c << ' ' << ss << '\n';
    }
    std::cout << fockspace.size() << '\n';

}
