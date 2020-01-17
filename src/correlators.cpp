#include "correlators.hpp"

#include "operator.hpp"


double computeCorrelatorNormalisation(Spectrum const &spectrum)
{
    double normalisation = 0.0;
    for (std::size_t i = 0; i < spectrum.size(); ++i) {
        normalisation += std::exp(-beta * spectrum.energies[i]);
    }
    return normalisation;
}


Correlators computeCorrelators(Spectrum const &spectrum)
{
    double const normalisation = computeCorrelatorNormalisation(spectrum);

    std::vector<DSparseMatrix> annihilatorElements;
    annihilatorElements.reserve(NSITES);
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
                /*
                 * This implements
                 *   Tr(B * A_i * C * A_j^T)
                 *   B = exp((tau-beta) * E)
                 *   C = exp(-tau * E)
                 */
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
