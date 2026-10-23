from root import find_roots
import numpy as np
import matplotlib.pyplot as plt
from eq73 import eq_73
from constants import *

def growth_to_frequency_contour(
    dispersion,
    kx_vals,
    ky_vals,
    run=0,
    save_path=None,
):
    """
    Plot Im(omega) / Re(omega) over the wave-vector grid.

    Uses the accepted root with the largest imaginary part at each grid point.
    The ratio is undefined where Re(omega) == 0, so those points are left blank.
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

    # Compute Im(omega) / Re(omega), masking zero denominators.
    real_omega = omega_grid.real
    imag_omega = omega_grid.imag
    ratio = np.full(shape, np.nan, dtype=float)

    valid = (
        np.isfinite(real_omega)
        & np.isfinite(imag_omega)
        & (real_omega != 0.0)
    )
    ratio[valid] = imag_omega[valid] / real_omega[valid]
    ratio = np.ma.masked_invalid(ratio)

    fig, ax = plt.subplots(figsize=(8, 6))

    if ratio.count() == 0:
        ax.text(
            0.5,
            0.5,
            "No valid roots found",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    else:
        mesh = ax.pcolormesh(
            kx_vals,
            ky_vals,
            ratio,
            shading="auto",
            cmap="RdBu_r",
        )

        fig.colorbar(
            mesh,
            ax=ax,
            label=r"$\mathrm{Im}(\omega)/\mathrm{Re}(\omega)$",
        )

    ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
    ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    ax.set_title(
        rf"$\mathrm{{Im}}(\omega)/\mathrm{{Re}}(\omega)$"
        f" — altitude {alt_eff[run]:g} km"
    )

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    plt.show()

    return omega_grid

kx_vals = np.linspace(-100, 100, 100)
ky_vals = np.linspace(-100, 100, 100)

growth_to_frequency_contour(
    eq_73,
    kx_vals,
    ky_vals,
    run=5,
    save_path="growth_to_frequency_104km.png",
)