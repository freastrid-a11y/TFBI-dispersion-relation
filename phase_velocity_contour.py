from root import find_roots
import numpy as np
import matplotlib.pyplot as plt
from eq73 import eq_73
from constants import *

def phase_velocity_contour(
    dispersion,
    kx_vals,
    ky_vals,
    run=0,
    save_path=None,
):
    """
    Plot Re(omega) / |k| over the wave-vector grid.

    Uses the accepted root with the largest imaginary part at each grid point,
    as contour_plot() does. The k = 0 point is left blank because the quantity
    is undefined there.
    """
    kx_vals = np.asarray(kx_vals, dtype=float)
    ky_vals = np.asarray(ky_vals, dtype=float)

    shape = (len(ky_vals), len(kx_vals))
    omega_grid = np.full(shape, np.nan + 1j * np.nan, dtype=complex)

    for iy, ky in enumerate(ky_vals):
        for ix, kx in enumerate(kx_vals):
            continuation_seeds = []

            # Seed from previously solved neighboring points.
            if ix > 0 and np.isfinite(omega_grid[iy, ix - 1]):
                continuation_seeds.append(omega_grid[iy, ix - 1])

            if iy > 0 and np.isfinite(omega_grid[iy - 1, ix]):
                continuation_seeds.append(omega_grid[iy - 1, ix])

            roots = find_roots(
                dispersion,
                kx,
                ky,
                run,
                guesses=continuation_seeds,
                include_default_guesses=True,
            )

            if roots:
                omega_grid[iy, ix] = roots[0]

        print(f"Altitude {alt_eff[run]:g} km: row {iy + 1}/{len(ky_vals)}")

    # Compute Re(omega) / sqrt(kx^2 + ky^2).
    KX, KY = np.meshgrid(kx_vals, ky_vals)
    k_mag = np.hypot(KX, KY)

    phase_velocity = np.full(shape, np.nan, dtype=float)
    valid = (
        np.isfinite(omega_grid.real)
        & np.isfinite(k_mag)
        & (k_mag > 0.0)
    )
    phase_velocity[valid] = omega_grid.real[valid] / k_mag[valid]
    phase_velocity = np.ma.masked_invalid(phase_velocity)

    fig, ax = plt.subplots(figsize=(8, 6))

    if phase_velocity.count() == 0:
        ax.text(
            0.5,
            0.5,
            "No accepted roots found",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    else:
        mesh = ax.pcolormesh(
            kx_vals,
            ky_vals,
            phase_velocity,
            shading="auto",
            cmap="RdBu_r",
        )

        fig.colorbar(
            mesh,
            ax=ax,
            label=r"$\mathrm{Re}(\omega)/\sqrt{k_x^2+k_y^2}$ [m s$^{-1}$]",
        )

    ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
    ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    ax.set_title(f"Phase velocity — altitude {alt_eff[run]:g} km")

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    plt.show()

    return omega_grid

kx_vals = np.linspace(-100, 100, 100)
ky_vals = np.linspace(-100, 100, 100)

phase_velocity_contour(
    eq_73,
    kx_vals,
    ky_vals,
    run=4,
    save_path="phase_velocity_110km.png",
)