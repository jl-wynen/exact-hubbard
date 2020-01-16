#include <iostream>
#include <fstream>
#include <optional>
#include <filesystem>
#include <chrono>

#include "operator.hpp"
#include "state.hpp"
#include "spectrum.hpp"

namespace fs = std::filesystem;


void saveSpectrum(fs::path const &fname, Spectrum const &spectrum)
{
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


    [[nodiscard]] static std::size_t totalIndex(std::size_t const i,
                                                std::size_t const j,
                                                std::size_t const t) noexcept
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
DSparseMatrix toEigenspaceMatrix(DMatrix const &matrix, Spectrum const &spectrum)
{
    DSparseMatrix res(matrix.rows(), matrix.columns());

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
            if (blaze::abs(elem) > 1e-8) {
                res(alpha, gamma) = elem;
            }
        }
    }

    return res;
}


Correlators computeCorrelators(Spectrum const &spectrum)
{
    double const normalisation = computeCorrelatorNormalisation(spectrum);

    std::vector<DSparseMatrix> annihilatorElements;
    for (std::size_t i = 0; i < NSITES; ++i) {
        DSparseMatrix const elementsFockspace = toMatrix(ParticleAnnihilator{i},
                                                   spectrum.basis);
        annihilatorElements.emplace_back(toEigenspaceMatrix(elementsFockspace, spectrum));
    }

    Correlators corrs;
    for (std::size_t i = 0; i < NSITES; ++i) {
        auto const &Ai = annihilatorElements[i];
        for (std::size_t j = 0; j < NSITES; ++j) {
            auto const &Aj = annihilatorElements[j];

            for (std::size_t t = 0; t < NT; ++t) {
                double const tau = beta / static_cast<double>(NT - 1) * static_cast<double>(t);
                corrs(i, j, t) = trace((expand(exp((tau-beta)*spectrum.energies), spectrum.size())
                                        % Ai)
                                       * (expand(exp(-tau*spectrum.energies), spectrum.size())
                                          % trans(Aj)))
                                 / normalisation;
            }
        }
    }

    return corrs;
}


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
    correlators.save("../correlators.dat");
}
