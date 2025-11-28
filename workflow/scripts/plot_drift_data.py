import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.style.use("styles/default.mplstyles")


def plot_npz_file(input_filename, output_prefix):
    npz = np.load(input_filename)
    gridx = npz["gridx"]
    gridy = npz["gridy"]
    score = npz["score"]
    logP = npz["logP"]
    P = npz["pdf"]

    U, V = score[..., 0], score[..., 1]
    x_lims = (gridx.min(), gridx.max())
    y_lims = (gridy.min(), gridy.max())

    if "EBM" in input_filename:
        # 1. potential surface
        fig = plt.figure(figsize=(9, 6), dpi=100)
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(
            gridx,
            gridy,
            np.clip(logP, 0, 10),
            cmap="Wistia",
            linewidth=0.4,
            edgecolor="k",
        )
        ax.set_xlim(*x_lims)
        ax.set_ylim(*y_lims)
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.view_init(elev=15, azim=20)
        ax.tick_params(axis="both", which="major", labelsize=10, direction="in")
        ax.grid(False)
        ax.set_box_aspect([1.2, 1, 0.4], zoom=1.2)  # Better shape
        plt.savefig(f"{output_prefix}_potential.pdf", bbox_inches="tight")
        plt.close()

        # 2. distribution
        fig = plt.figure(figsize=(9, 6), dpi=100)
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(gridx, gridy, P, cmap="Wistia", linewidth=0.4, edgecolor="k")
        ax.set_xlim(*x_lims)
        ax.set_ylim(*y_lims)
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.view_init(elev=15, azim=20)
        ax.tick_params(axis="both", which="major", labelsize=10, direction="in")
        ax.grid(False)
        ax.set_box_aspect([1.2, 1, 0.4], zoom=1.2)  # Better shape
        plt.savefig(f"{output_prefix}_distribution.pdf", bbox_inches="tight")
        plt.close()

    # 3. drift streamplot
    fig, ax = plt.subplots(figsize=(9, 6), dpi=100)
    ax.streamplot(
        gridx[0, :],
        gridy[:, 0],
        U,
        V,
        color="b",
        linewidth=1,
        arrowstyle="->",
        density=2,
    )
    ax.set_xlim(*x_lims)
    ax.set_ylim(*y_lims)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    plt.savefig(f"{output_prefix}_streamplot.pdf", bbox_inches="tight")
    plt.close()

    # Thimbles only for QMODEL
    if "QMODEL" in input_filename:
        x = gridx[0, :]
        D3 = 1 + 18 * x * x * (1 - 2 * x * x)
        D2 = 1 + 12 * x * x * (1 + x * x)
        D1 = (D3 + np.sqrt(D3**2 - D2**3 + 0j)) ** (1 / 3)

        phase = np.exp(-1j * np.pi / 3)
        stable_thimble = x + 1j * (-1 + phase * (D2 / D1) + D1 / phase) / (6 * x + 1e-6)

        phase = np.exp(1j * np.pi / 3)
        unstable_thimble = x + 1j * (-1 + phase * (D2 / D1) + D1 / phase) / (
            6 * x + 1e-6
        )

        fig, ax = plt.subplots(figsize=(9, 6), dpi=100)
        ax.streamplot(
            gridx[0, :],
            gridy[:, 0],
            U,
            V,
            density=2,
            color="b",
            linewidth=1,
            arrowstyle="->",
        )
        ax.plot(stable_thimble.real, stable_thimble.imag, "k", lw=3)
        ax.plot(unstable_thimble.real, unstable_thimble.imag, "r--", lw=2)
        ax.set_xlim(*x_lims)
        ax.set_ylim(*y_lims)
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        plt.savefig(f"{output_prefix}_streamplot_w_thimbles.pdf", bbox_inches="tight")
        plt.close()


if __name__ == "__main__":
    plot_npz_file(
        input_filename=snakemake.input[0], output_prefix=snakemake.params.prefix
    )
