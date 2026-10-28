"""Find the maximum growth_rate / frequency ratio in a growth-rate CSV and plot it."""

import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from constants import alt_eff, nu_i

FILES = [
    "data/104.csv",
    "data/110.csv",
]


def load_ratio_grid(path):
    """Return (kx_vals, ky_vals, ratio_grid, omega_mag_grid) pivoted from the long-format CSV."""
    df = pd.read_csv(path)
    ratio = df["growth_rate_s^-1"] / df["frequency_rad_s^-1"]
    df["ratio"] = ratio.replace([np.inf, -np.inf], np.nan)
    df["omega_mag"] = np.hypot(df["growth_rate_s^-1"], df["frequency_rad_s^-1"])

    pivot_ratio = df.pivot(index="ky_m^-1", columns="kx_m^-1", values="ratio")
    pivot_omega = df.pivot(index="ky_m^-1", columns="kx_m^-1", values="omega_mag")
    return (
        pivot_ratio.columns.to_numpy(),
        pivot_ratio.index.to_numpy(),
        pivot_ratio.to_numpy(),
        pivot_omega.to_numpy(),
    )


def max_ratio(path):
    _, _, ratio, _ = load_ratio_grid(path)
    return np.nanmax(np.abs(ratio))


def nu_in_for_altitude(path):
    """Look up nu_i for the altitude (km) encoded in the filename, using the first 3-digit run found, e.g. '110_smallk.csv' -> alt_eff==110."""
    altitude = int(re.search(r"\d{3}", path.split("/")[-1]).group())
    return nu_i[np.where(alt_eff == altitude)[0][0]]


def plot_ratio_contour(path, save_path, band_frac=0.2):
    """Growth/frequency ratio contour, color scale fixed symmetric about this file's max |ratio|.

    band_frac sets how wide the accepted |omega|~nu_in resonance band is, as a
    fraction of nu_in (e.g. 0.2 -> accept |omega| within +/-20% of nu_in).
    """
    kx_vals, ky_vals, ratio, omega_mag = load_ratio_grid(path)

    nu_in = nu_in_for_altitude(path)
    lower = nu_in * (1.0 - band_frac)
    upper = nu_in * (1.0 + band_frac)

    in_band = (
        np.isfinite(omega_mag)
        & (omega_mag >= lower)
        & (omega_mag <= upper)
        & np.isfinite(ratio)
    )

    if not np.any(in_band):
        print(f"No sampled points in the transition band: {path}")
        return

    # Keep only ratios at sampled points inside the band.
    band_ratio = np.ma.masked_where(~in_band, ratio)

    abs_values = np.abs(ratio[in_band])

    print(
        f"Transition-band points: {in_band.sum()}\n"
        f"Median |Im/Re|: {np.median(abs_values):.4g}\n"
        f"Maximum |Im/Re|: {np.max(abs_values):.4g}"
    )

    # A percentile limit prevents a few outliers from hiding the structure.
    vmax = max(float(np.percentile(abs_values, 95)), 1e-12)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(8, 6))
    mesh = ax.pcolormesh(
        kx_vals,
        ky_vals,
        band_ratio,
        shading="auto",
        cmap="RdBu_r",
        norm=norm,
    )

    fig.colorbar(
        mesh,
        ax=ax,
        extend="both",
        label=r"$\mathrm{Im}(\omega)/\mathrm{Re}(\omega)$",
    )

    # Mark the resonance band where |omega| is within band_frac of nu_in (ion-neutral collision frequency).
    nu_in = nu_in_for_altitude(path)
    lower, upper = nu_in * (1.0 - band_frac), nu_in * (1.0 + band_frac)
    omega_mag = np.ma.masked_invalid(omega_mag)
    found = bool(omega_mag.count() and omega_mag.min() < upper and omega_mag.max() > lower)
    if found:
        ax.contourf(
            kx_vals, ky_vals, omega_mag, levels=[lower, upper],
            colors="none", hatches=["//"],
        )
        contour = ax.contour(kx_vals, ky_vals, omega_mag, levels=[lower, upper], colors="black", linewidths=1.0)
        ax.clabel(contour, fmt={lower: r"$-$", upper: r"$+$"}, fontsize=8, inline=True)
        print(
            f"  |omega|~nu_in band found: nu_in={nu_in:.3g} (band=[{lower:.3g}, {upper:.3g}]), "
            f"|omega| range=[{omega_mag.min():.3g}, {omega_mag.max():.3g}]"
        )
    else:
        range_str = f"[{omega_mag.min():.3g}, {omega_mag.max():.3g}]" if omega_mag.count() else "empty"
        print(f"  no |omega|~nu_in band: nu_in band=[{lower:.3g}, {upper:.3g}] not within |omega| range={range_str}")

    ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
    ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    ax.set_title(f"Growth/frequency ratio (|max| = {vmax:.3g})")
    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {save_path}")


def load_field_grid(path, column="growth_rate_s^-1"):
    """Return (kx_vals, ky_vals, grid) pivoted from the long-format CSV for the given column."""
    df = pd.read_csv(path)
    pivot = df.pivot(index="ky_m^-1", columns="kx_m^-1", values=column)
    return pivot.columns.to_numpy(), pivot.index.to_numpy(), pivot.to_numpy()


def plot_field_difference(path_nonzero, path_zero, save_path, column="growth_rate_s^-1"):
    """Contour of path_nonzero minus path_zero (e.g. a_i0 != 0 vs a_i0 == 0), color scale symmetric about the max |difference|.

    Assumes both CSVs were sampled on the same kx/ky grid.
    """
    kx_vals, ky_vals, grid_nonzero = load_field_grid(path_nonzero, column)
    kx_vals_zero, ky_vals_zero, grid_zero = load_field_grid(path_zero, column)

    if not (np.array_equal(kx_vals, kx_vals_zero) and np.array_equal(ky_vals, ky_vals_zero)):
        raise ValueError("kx/ky grids differ between the two files; cannot take a pointwise difference.")

    diff = grid_nonzero - grid_zero
    vmax = max(float(np.nanmax(np.abs(diff))), 1e-12)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(8, 6))
    mesh = ax.pcolormesh(kx_vals, ky_vals, diff, shading="auto", cmap="RdBu_r", norm=norm)
    fig.colorbar(mesh, ax=ax, label=rf"$\Delta$ {column}")

    ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
    ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    ax.set_title(f"{column} difference (a_i0$\\neq$0 minus a_i0=0), max |diff| = {vmax:.3g}")
    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {save_path}")


if __name__ == "__main__":
    for path in FILES:
        print(path, max_ratio(path))
        save_path = "plots/" + path.split("/")[-1].replace(".csv", "_ratio.png")
        plot_ratio_contour(path, save_path)
