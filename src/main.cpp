#include <iostream>

#include "state.hpp"
#include "io.hpp"
#include "operator.hpp"


int main()
{
    auto fockspace = fockspaceBasis();

    for (auto s : fockspace) {
        std::cout << s << '\n';
    }
    std::cout << fockspace.size() << '\n';

}
