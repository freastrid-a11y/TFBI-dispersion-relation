# plot_growth_damping.py
#
# Loads the .npz produced by generate_growth_damping.py and plots it.
# Edit anything visual here and re-run without recomputing the roots.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


# ---------------- SETTINGS ----------------
RUN = 4

INPUT_FILE = f"dispersion_run_{RUN}_roots.npz"
# ------------------------------------------


def load_data():
    with np.load(INPUT_FILE) as data:
        return (
            data["kx"],
            data["ky"],
            data["omega"],
            float(data["altitude"]),
        )


def symmetric_norm(data):
    values = np.ma.masked_invalid(data).compressed()

    limit = (
        max(float(np.max(np.abs(values))), 1e-12)
        if values.size else 1.0
    )

    return TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)

def actual_range_norm(data):
    values = np.ma.masked_invalid(data).compressed()

    if values.size == 0:
        return TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)

    minimum = float(np.min(values))
    maximum = float(np.max(values))

    # TwoSlopeNorm requires vmin < 0 < vmax.
    return TwoSlopeNorm(
        vmin=min(minimum, -1e-12),
        vcenter=0.0,
        vmax=max(maximum, 1e-12),
    )


def add_zero_growth_boundary(ax, kx_values, ky_values, gamma):
    """Overlay the marginal-stability contour wherever it can be resolved."""
    gamma = np.ma.masked_invalid(gamma)

    if gamma.count() and gamma.min() < 0.0 < gamma.max():
        contour = ax.contour(
            kx_values,
            ky_values,
            gamma,
            levels=[0.0],
            colors="black",
            linewidths=1.4,
            corner_mask=False,
        )

        ax.clabel(
            contour,
            fmt={0.0: r"$\gamma=0$"},
            fontsize=8,
            inline=True,
        )


def plot_results(kx_values, ky_values, omega_grid, altitude):
    valid = np.isfinite(omega_grid)

    gamma = np.where(valid, omega_grid.imag, np.nan)
    omega_real = np.where(valid, omega_grid.real, np.nan)

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("lightgray")

    fig, axes = plt.subplots(
        1, 2,
        figsize=(13, 6),
        constrained_layout=True,
    )

    panels = [
        (
            gamma,
            actual_range_norm(gamma),
            r"$\mathrm{Im}(\omega)$",
            r"$\mathrm{Im}(\omega)$ [s$^{-1}$]",
        ),
        (
            omega_real,
            actual_range_norm(omega_real),
            r"$\mathrm{Re}(\omega)$",
            r"$\mathrm{Re}(\omega)$ [s$^{-1}$]",
        ),
    ]

    for ax, panel in zip(axes, panels):
        data, norm, title, colorbar_label = panel

        mesh = ax.pcolormesh(
            kx_values,
            ky_values,
            np.ma.masked_invalid(data),
            shading="auto",
            cmap=cmap,
            norm=norm,
        )

        cbar = fig.colorbar(
            mesh,
            ax=ax,
            label=colorbar_label,
        )

        # Four intervals on each side of zero.
        negative_ticks = np.linspace(norm.vmin, 0.0, 5)
        positive_ticks = np.linspace(0.0, norm.vmax, 5)[1:]

        ticks = np.concatenate((negative_ticks, positive_ticks))

        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f"{value:,.0f}" for value in ticks])

        add_zero_growth_boundary(ax, kx_values, ky_values, gamma)

        ax.set_title(title)
        ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
        ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
        ax.set_aspect("equal")

    fig.suptitle(
        f"Effective altitude {altitude:g} km\n",
        fontsize=14,
    )

    filename = f"dispersion_run_{RUN}_growth_damping.png"
    fig.savefig(filename, dpi=250)
    print(f"Saved figure: {filename}")

    plt.show()


def main():
    kx_values, ky_values, omega_grid, altitude = load_data()
    plot_results(kx_values, ky_values, omega_grid, altitude)


if __name__ == "__main__":
    main()