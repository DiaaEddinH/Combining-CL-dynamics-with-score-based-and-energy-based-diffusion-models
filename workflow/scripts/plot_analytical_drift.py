import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from typing import Optional, Tuple

plt.style.use("styles/default.mplstyles")


def plot_analytical_quartic_drift_streamplot(
    output_filename: str,
    mass: complex,
    coupling: complex,
    no_of_points: int = 101,
    x_lims: Tuple[float] = (-3, 3),
    y_lims: Tuple[float] = (-3, 3),
):
    x = np.linspace(*x_lims, no_of_points)
    y = np.linspace(*y_lims, no_of_points)
    gridx, gridy = np.meshgrid(x, y)

    z = gridx + 1j * gridy

    F = -mass * z - coupling * z**3

    U, V = F.real, F.imag

    fig, ax = plt.subplots(1, 1, dpi=100, figsize=(9, 6))

    ax.streamplot(
        gridx, gridy, U, V, linewidth=1, color="b", arrowstyle="->", density=2
    )
    ax.set_xlim(*x_lims)
    ax.set_ylim(*y_lims)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches="tight")
    plt.close(fig)


def plot_analytical_gaussian_drift_streamplot(
    output_filename: str,
    mu: tuple,
    sigma: float,
    no_of_points: int = 101,
    x_lims: Tuple[float] = (-3, 3),
    y_lims: Tuple[float] = (-3, 3),
):
    mu = np.array(mu)
    x = np.linspace(*x_lims, no_of_points)
    y = np.linspace(*y_lims, no_of_points)
    gridx, gridy = np.meshgrid(x, y)

    X = np.stack([gridx, gridy], axis=-1)

    F = (
        -(
            X
            - mu[None, None]
            * np.tanh((mu[0] * X[..., 0] + mu[1] * X[..., 1])[..., None] / sigma**2)
        )
        / sigma**2
    )
    U, V = F[..., 0], F[..., 1]

    fig, ax = plt.subplots(1, 1, dpi=100, figsize=(9, 6))

    ax.streamplot(
        gridx, gridy, U, V, linewidth=1, color="b", arrowstyle="->", density=2
    )
    ax.set_xlim(*x_lims)
    ax.set_ylim(*y_lims)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    plot_analytical_quartic_drift_streamplot(
        output_filename=snakemake.output.pdf_cl,
        mass=complex(snakemake.params.mass),
        coupling=snakemake.params.coupling,
        x_lims=snakemake.params.x_lims,
        y_lims=snakemake.params.y_lims,
    )

    plot_analytical_gaussian_drift_streamplot(
        output_filename=snakemake.output.pdf_gaussian,
        mu=snakemake.params.mu,
        sigma=snakemake.params.sigma,
        x_lims=(-3, 3),
        y_lims=(-3, 3),
    )
