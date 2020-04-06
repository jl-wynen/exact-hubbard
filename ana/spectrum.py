from collections import Counter

import numpy as np
import matplotlib.pyplot as plt


def collect_degenerates(spectrum):
    q_to_e = dict()
    for charge, energy in spectrum:
        if charge in q_to_e:
            q_to_e[charge].append(energy)
        else:
            q_to_e[charge] = [energy]

    q_to_e = {charge: Counter(energies) for charge, energies in q_to_e.items()}
    charges, energies = zip(*sorted(q_to_e.items(), key=lambda t: t[0]))
    return charges, energies


def main():
    spectrum = np.loadtxt("spectrum.dat")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("Q")
    ax.set_ylabel(r"$(E_Q^\alpha - E_\Omega) / \kappa$")

    spectrum[:, 1] = np.round(spectrum[:, 1] - np.min(spectrum[:, 1]), 5)
    for charge, energies in zip(*collect_degenerates(spectrum)):
        for energy, count in energies.items():
            ax.plot((charge-0.33, charge+0.33), [energy]*2,
                    linewidth=2, c="C0")
            ax.text(charge, energy+0.1, f"({count})", horizontalalignment="center")

    fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
