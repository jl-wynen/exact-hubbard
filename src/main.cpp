#include <iostream>
#include <fstream>
#include <optional>
#include <filesystem>

#include "operator.hpp"
#include "state.hpp"
#include "spectrum.hpp"

namespace fs = std::filesystem;


/**
 * Compute and store the non-interacting spectrum.
 * That is, compute the simultaneous eigen system of
 * the Hamiltonian and charge operator.
 * The result is written to a file named `fname`.
 */
void computeNonInteractingSpectrum(fs::path const &fname="../spectrum.dat")
{
    auto const fockspace = fockspaceBasis();
    auto const spectrum = Spectrum::compute(fockspace);
    std::ofstream ofs(fname);
    ofs << "#  Q  E\n";
    for (std::size_t i = 0; i < spectrum.size(); ++i) {
        ofs << spectrum.charges[i] << ' ' << spectrum.energies[i] << '\n';
    }
}


double computeCorrelatorNormalisation(Spectrum const &spectrum)
{
    double normalisation = 0.0;
    for (std::size_t i = 0; i < spectrum.size(); ++i) {
        normalisation += std::exp(-beta * spectrum.energies[i]);
    }
    return normalisation;
}


class Correlators
{
    std::vector<double> data_;

public:
    explicit Correlators()
            : data_(NSITES*NSITES*NT)
    { }


    [[nodiscard]] std::size_t totalIndex(std::size_t const i,
                                         std::size_t const j,
                                         std::size_t const t) const noexcept
    {
        assert(i < NSITES);
        assert(j < NSITES);
        assert(t < NT);
        return (i*NSITES + j)*NT + t;
    }


    double operator()(std::size_t const i,
                      std::size_t const j,
                      std::size_t const t) const noexcept
    {
        return data_[totalIndex(i, j, t)];
    }


    double &operator()(std::size_t const i,
                       std::size_t const j,
                       std::size_t const t) noexcept
    {
        return data_[totalIndex(i, j, t)];
    }


    void save(fs::path const &fname) const
    {
        std::ofstream ofs{fname};
        ofs << "#~ correlator\n#  nx  nt\n"
            << NSITES << ' ' << NT
            << "\n#  U  kappa  beta\n"
            << U << ' ' << kappa << ' ' << beta
            << "\n#  data\n";
        for (auto const x : data_) {
            ofs << x << ' ';
        }
    }
};


/**
 * Turn matrix elements for an operator in basis `spectrum.basis`
 * into matrix elements in the basis of eigenvectors.
 * \param matrix Matrix elements in `spectrum.basis`.
 * \param spectrum Provides basis and eigenstates.
 * \return Matrix elements in eigenbasis.
 */
DMatrix toEigenspaceMatrix(DMatrix const &matrix, Spectrum const &spectrum)
{
    DMatrix res(matrix.rows(), matrix.columns());

    for (std::size_t alpha = 0; alpha < spectrum.size(); ++alpha) {
        for (std::size_t gamma = 0; gamma < spectrum.size(); ++gamma) {
            double elem = 0.0;
            for (std::size_t x = 0; x < spectrum.eigenStateIdxs[alpha].size(); ++x) {
                for (std::size_t y = 0; y < spectrum.eigenStateIdxs[gamma].size(); ++y) {
                    elem += spectrum.eigenStateCoeffs[alpha][x]
                            * spectrum.eigenStateCoeffs[gamma][y]
                            * matrix(spectrum.eigenStateIdxs[alpha][x],
                                     spectrum.eigenStateIdxs[gamma][y]);
                }
            }
            res(alpha, gamma) = elem;
        }
    }

    return res;
}


void computeInteractingCorrelators(fs::path const &fname="../correlators.dat")
{
    auto const spectrum = Spectrum::compute(fockspaceBasis());
    double const normalisation = computeCorrelatorNormalisation(spectrum);

    std::vector<DMatrix> annihilatorElements;
    for (std::size_t i = 0; i < NSITES; ++i) {
        DMatrix const elementsFockspace = toMatrix(ParticleAnnihilator{i},
                                                   spectrum.basis);
        annihilatorElements.emplace_back(toEigenspaceMatrix(elementsFockspace, spectrum));
    }

    Correlators corrs;
    for (std::size_t i = 0; i < NSITES; ++i) {
        DMatrix const &Ai = annihilatorElements[i];
        for (std::size_t j = 0; j < NSITES; ++j) {
            DMatrix const &Aj = annihilatorElements[j];

            for (std::size_t t = 0; t < NT; ++t) {
                double const tau = beta / static_cast<double>(NT - 1) * static_cast<double>(t);
                double accum = 0.0;
                for (std::size_t alpha = 0; alpha < Ai.rows(); ++alpha) {
                    for (std::size_t gamma = 0; gamma < Ai.columns(); ++gamma) {
                        accum += std::exp((tau-beta)*spectrum.energies[alpha] - tau*spectrum.energies[gamma])
                                 * Ai(alpha, gamma) * Aj(alpha, gamma);
                    }
                }
                corrs(i, j, t) = accum / normalisation;
            }
        }
    }

    corrs.save(fname);
}


int main()
{
   // computeNonInteractingSpectrum();
   computeInteractingCorrelators();
}
