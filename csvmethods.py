"""Plot the point-by-point difference between two saved grid CSVs in data/."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


def load_field_grid(path, column="growth_rate_s^-1"):
    """Read one column of a grid CSV; return kx, ky and the 2D grid."""
    df = pd.read_csv(path)
    pivot = df.pivot(index="ky_m^-1", columns="kx_m^-1", values=column)
    return pivot.columns.to_numpy(), pivot.index.to_numpy(), pivot.to_numpy()


def plot_field_difference(path_nonzero, path_zero, save_path, column="growth_rate_s^-1", percentile=95):
    """
    Plot path_nonzero minus path_zero, e.g. (a_i0 on) minus (a_i0 off).

    Both CSVs must use the same k grid. The colour scale is symmetric
    about 0 and set by the given percentile of |difference|, so a few
    extreme points don't wash out the plot; the colour-bar arrows mark
    values beyond that.
    """
    kx_vals, ky_vals, grid_nonzero = load_field_grid(path_nonzero, column)
    kx_vals_zero, ky_vals_zero, grid_zero = load_field_grid(path_zero, column)

    if not (np.array_equal(kx_vals, kx_vals_zero) and np.array_equal(ky_vals, ky_vals_zero)):
        raise ValueError("kx/ky grids differ between the two files; cannot take a pointwise difference.")

    diff = grid_nonzero - grid_zero
    abs_diff = np.abs(diff)
    finite_abs_diff = abs_diff[np.isfinite(abs_diff)]
    if finite_abs_diff.size == 0:
        raise ValueError("No finite difference values to plot.")

    vmax = max(float(np.nanpercentile(finite_abs_diff, percentile)), 1e-12)
    actual_max = float(np.nanmax(finite_abs_diff))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(8, 6))
    mesh = ax.pcolormesh(kx_vals, ky_vals, diff, shading="auto", cmap="RdBu_r", norm=norm)
    fig.colorbar(mesh, ax=ax, extend="both", label=rf"$\Delta$ {column}")

    ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
    ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    ax.set_title(
        f"{column} difference (a_i0$\\neq$0 minus a_i0=0),",
    )
    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {save_path}")


if __name__ == "__main__":
    plot_field_difference("data/104.csv", "data/104_ai0zero.csv", "plots/difference_104")
    plot_field_difference("data/110.csv", "data/110_ai0zero.csv", "plots/difference_110")
