#include <iostream>
#include <iomanip>
#include <fstream>

#include "state.hpp"
#include "io.hpp"
#include "operator.hpp"
#include "matrix.hpp"


template <typename OP>
void computeSpectrum(OP const &hamiltonian,
                     SumState const &basis)
{
    auto const mat = toMatrix(hamiltonian, basis);
    std::ofstream ofs("../hamiltonian.mat");
    ofs << mat << '\n';
};


int main()
{
    auto const fockspace = fockspaceBasis();

    SumOperator hamiltonian{ParticleHop{},
                            HoleHop{},
                            SquaredNumberOperator{}};

    computeSpectrum(hamiltonian, fockspace);
}
