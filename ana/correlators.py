import numpy as np
import matplotlib.pyplot as plt


def linestyle(i):
    linestyles = ["-", "--", "-.", ":"]
    return linestyles[i]


def get_irreps(kappa):
    """
    Compute the lattice irreps.
    """

    # Square
    # hopping = np.array([[0, 1, 0, 1],
    #                     [1, 0, 1, 0],
    #                     [0, 1, 0, 1],
    #                     [1, 0, 1, 0]])

    # # Triangle
    # hopping = np.array([[0, 1, 1],
    #                     [1, 0, 1],
    #                     [1, 1, 0]])

    #  tetrahedron
    hopping = np.array([[0, 1, 1, 1],
                        [1, 0, 1, 1],
                        [1, 1, 0, 1],
                        [1, 1, 1, 0]])

    return np.linalg.eigh(hopping * kappa)[1]


def project_to_irreps(corrs, params):
    """
    Project correlators from position space to the irrep basis.
    """
    irreps = get_irreps(params["kappa"])
    return np.einsum("ij,jkt,kl->ilt", irreps.T, corrs, irreps)


def load_correlators(fname):
    """
    Load correlators and meta data stored in a file.
    """

    with open(fname, "r") as f:
        assert f.readline() == "#~ correlator\n"
        assert f.readline() == "#  nx  nt\n"
        nx, nt = map(int, f.readline().split(" "))
        assert f.readline() == "#  U  kappa  beta\n"
        U, kappa, beta = map(float, f.readline().split(" "))

    corrs = np.loadtxt(fname, skiprows=5).reshape(nx, nx, nt)

    return corrs, dict(U=U, kappa=kappa, beta=beta)


def plot_all_in_one(corrs, params):
    """
    Plot all correlators in one plot.
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(r"$\kappa \tau$")
    ax.set_ylabel(r"$C(\tau)$")

    x = np.linspace(0, params["beta"], corrs.shape[2], endpoint=True) * params["kappa"]
    for i, j in np.ndindex(corrs.shape[:2]):
        ax.plot(x, corrs[i, j], c=f"C{i}", ls=linestyle(j))

    fig.tight_layout()


def plot_grid(corrs, params):
    """
    Plot a grid of all correlators.
    """

    fig = plt.figure(figsize=(11, 10))
    fig.suptitle(rf"$U/\kappa = {params['U']/params['kappa']} \qquad \kappa \beta = {params['kappa']*params['beta']}$")

    x = np.linspace(0, params["beta"], corrs.shape[2], endpoint=True) * params["kappa"]
    for i, j in np.ndindex(corrs.shape[:2]):
        ax = fig.add_subplot(corrs.shape[0], corrs.shape[1], i*corrs.shape[1] + j + 1)
        ax.set_xlabel(r"$\kappa \tau$")
        ax.set_ylabel(rf"$C_{{{i},{j}}}(\tau)$")
        ax.plot(x, corrs[i, j])
        ax.set_yscale("log")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])


def main():
    corrs, params = load_correlators("../correlators.dat")
    corrs = project_to_irreps(corrs, params)

    plot_grid(corrs, params)
    plt.show()


if __name__ == '__main__':
    main()
