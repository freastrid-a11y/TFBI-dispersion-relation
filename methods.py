"""
Plotting methods that consume an already-solved (kx, ky, growth, frequency)
grid, as produced/cached by grid_cache.get_grid(). Each method has the
signature (kx_vals, ky_vals, growth, frequency, altitude, save_path) and
saves a figure to save_path.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


def _actual_range_norm(data):
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


def _add_zero_growth_boundary(ax, kx_vals, ky_vals, growth):
    """Overlay the marginal-stability contour wherever it can be resolved."""
    growth = np.ma.masked_invalid(growth)

    if growth.count() and growth.min() < 0.0 < growth.max():
        contour = ax.contour(
            kx_vals,
            ky_vals,
            growth,
            levels=[0.0],
            colors="black",
            linewidths=1.4,
            corner_mask=False,
        )

        ax.clabel(contour, fmt={0.0: r"$\gamma=0$"}, fontsize=8, inline=True)


def growth_damping(kx_vals, ky_vals, growth, frequency, altitude, save_path):
    """Two-panel Im(omega) / Re(omega) contour."""
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("lightgray")

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)

    panels = [
        (growth, r"$\mathrm{Im}(\omega)$", r"$\mathrm{Im}(\omega)$ [s$^{-1}$]"),
        (frequency, r"$\mathrm{Re}(\omega)$", r"$\mathrm{Re}(\omega)$ [s$^{-1}$]"),
    ]

    for ax, (data, title, colorbar_label) in zip(axes, panels):
        norm = _actual_range_norm(data)

        mesh = ax.pcolormesh(
            kx_vals,
            ky_vals,
            np.ma.masked_invalid(data),
            shading="auto",
            cmap=cmap,
            norm=norm,
        )

        cbar = fig.colorbar(mesh, ax=ax, label=colorbar_label)

        # Four intervals on each side of zero.
        negative_ticks = np.linspace(norm.vmin, 0.0, 5)
        positive_ticks = np.linspace(0.0, norm.vmax, 5)[1:]
        ticks = np.concatenate((negative_ticks, positive_ticks))

        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f"{value:,.0f}" for value in ticks])

        _add_zero_growth_boundary(ax, kx_vals, ky_vals, growth)

        ax.set_title(title)
        ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
        ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
        ax.set_aspect("equal")

    fig.suptitle(f"Effective altitude {altitude:g} km\n", fontsize=14)
    fig.savefig(save_path, dpi=250)
    plt.close(fig)
    print(f"Saved figure: {save_path}")


def growth_contour(kx_vals, ky_vals, growth, frequency, altitude, save_path, positive_only=True):
    """Single-panel growth-rate (Im(omega)) contour."""
    data = growth.copy()

    if positive_only:
        data[data <= 0.0] = np.nan

    data = np.ma.masked_invalid(data)

    fig, ax = plt.subplots(figsize=(8, 6))

    if data.count() == 0:
        ax.text(
            0.5,
            0.5,
            "No accepted unstable roots found" if positive_only else "No accepted roots found",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    else:
        mesh = ax.pcolormesh(kx_vals, ky_vals, data, shading="auto", cmap="viridis")
        fig.colorbar(mesh, ax=ax, label=r"Growth rate $\mathrm{Im}(\omega)$ [s$^{-1}$]")

    ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
    ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    ax.set_title(f"Growth rate — altitude {altitude:g} km")
    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {save_path}")


def growth_frequency_ratio(kx_vals, ky_vals, growth, frequency, altitude, save_path):
    """Im(omega) / Re(omega) contour, undefined where Re(omega) == 0."""
    ratio = np.full(growth.shape, np.nan, dtype=float)

    valid = np.isfinite(growth) & np.isfinite(frequency) & (frequency != 0.0)
    ratio[valid] = growth[valid] / frequency[valid]
    ratio = np.ma.masked_invalid(ratio)

    fig, ax = plt.subplots(figsize=(8, 6))

    if ratio.count() == 0:
        ax.text(0.5, 0.5, "No valid roots found", ha="center", va="center", transform=ax.transAxes)
    else:
        mesh = ax.pcolormesh(kx_vals, ky_vals, ratio, shading="auto", cmap="RdBu_r")
        fig.colorbar(mesh, ax=ax, label=r"$\mathrm{Im}(\omega)/\mathrm{Re}(\omega)$")

    ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
    ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    ax.set_title(
        rf"$\mathrm{{Im}}(\omega)/\mathrm{{Re}}(\omega)$ — altitude {altitude:g} km"
    )
    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {save_path}")


def phase_velocity(kx_vals, ky_vals, growth, frequency, altitude, save_path):
    """Re(omega) / |k| contour, undefined at k = 0."""
    KX, KY = np.meshgrid(kx_vals, ky_vals)
    k_mag = np.hypot(KX, KY)

    result = np.full(growth.shape, np.nan, dtype=float)
    valid = np.isfinite(frequency) & (k_mag > 0.0)
    result[valid] = frequency[valid] / k_mag[valid]
    result = np.ma.masked_invalid(result)

    fig, ax = plt.subplots(figsize=(8, 6))

    if result.count() == 0:
        ax.text(0.5, 0.5, "No accepted roots found", ha="center", va="center", transform=ax.transAxes)
    else:
        mesh = ax.pcolormesh(kx_vals, ky_vals, result, shading="auto", cmap="RdBu_r")
        fig.colorbar(mesh, ax=ax, label=r"$\mathrm{Re}(\omega)/\sqrt{k_x^2+k_y^2}$ [m s$^{-1}$]")

    ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
    ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    ax.set_title(f"Phase velocity — altitude {altitude:g} km")
    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {save_path}")


METHODS = {
    "growth_damping": growth_damping,
    "growth_contour": growth_contour,
    "growth_frequency_ratio": growth_frequency_ratio,
    "phase_velocity": phase_velocity,
}
