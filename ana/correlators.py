from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def linestyle(i):
    linestyles = ["-", "--", "-.", ":"]
    return linestyles[i]


def get_irreps(kappa):
    hopping = np.array([[0, 1, 0, 1],
                        [1, 0, 1, 0],
                        [0, 1, 0, 1],
                        [1, 0, 1, 0]]) * kappa
    return np.linalg.eigh(hopping)[1]


def project_to_irreps(corrs, params):
    irreps = get_irreps(params["kappa"])
    res = np.empty_like(corrs)
    for t in range(corrs.shape[2]):
        res[:, :, t] = irreps @ corrs[:, :, t] @ irreps.T
    return res


def load_correlators(fname):
    with open(fname, "r") as f:
        assert f.readline() == "#~ correlator\n"
        assert f.readline() == "#  nx  nt\n"
        nx, nt = map(int, f.readline().split(" "))
        assert f.readline() == "#  U  kappa  beta\n"
        U, kappa, beta = map(float, f.readline().split(" "))

    corrs = np.loadtxt(fname, skiprows=5).reshape(nx, nx, nt)

    return corrs, dict(U=U, kappa=kappa, beta=beta)


def load_toms_correlators(nx, params):
    fname = Path("/home/jl/Uni/PhD/cns/ergodicity/plotting/figures/exactCorrelators") \
            / "4site_square.exact.U4.beta4.dat"
            # / f"{nx}site.exact.U{int(params['U'])}.beta{int(params['beta'])}.dat"
    data = np.loadtxt(fname, skiprows=1).T
    return data[0], data[1:]


def plot_all_in_one(corrs, params):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(r"$\kappa \tau$")
    ax.set_ylabel(r"$C(\tau)$")

    x = np.linspace(0, params["beta"], corrs.shape[2], endpoint=True) * params["kappa"]
    for i, j in np.ndindex(corrs.shape[:2]):
        ax.plot(x, corrs[i, j], c=f"C{i}", ls=linestyle(j))

    fig.tight_layout()


def plot_grid(corrs, params, tomtau, tomcorrs):
    fig = plt.figure()
    fig.suptitle(rf"$U/\kappa = {params['U']/params['kappa']} \qquad \kappa \beta = {params['kappa']*params['beta']}$")

    x = np.linspace(0, params["beta"], corrs.shape[2], endpoint=True) * params["kappa"]
    for i, j in np.ndindex(corrs.shape[:2]):
        ax = fig.add_subplot(corrs.shape[0], corrs.shape[1], i*corrs.shape[1] + j + 1)
        ax.set_xlabel(r"$\kappa \tau$")
        ax.set_ylabel(rf"$C_{{{i},{j}}}(\tau)$")
        ax.plot(x, corrs[i, j])
        if i == j:
            ax.plot(tomtau, tomcorrs[i], ls=":")
        ax.set_yscale("log")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])


def main():
    corrs, params = load_correlators("../correlators.dat")
    corrs = project_to_irreps(corrs, params)

    tomtau, tomcorrs, = load_toms_correlators(corrs.shape[0], params)

    plot_grid(corrs, params, tomtau, tomcorrs)
    plt.show()


if __name__ == '__main__':
    main()
