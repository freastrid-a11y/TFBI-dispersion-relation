# generate_growth_damping.py
#
# Solves the Eq. 73 dispersion relation on a kx/ky grid and saves the
# result to an .npz file. Run this whenever RUN, POINTS, or the k-ranges
# change. Use plot_growth_damping.py to (re-)plot the saved data.

import numpy as np

from eq73 import eq_73
from root import find_roots, alt_eff


# ---------------- SETTINGS ----------------
RUN = 5
POINTS = 100

KX_VALUES = np.linspace(-100.0, 100.0, POINTS)
KY_VALUES = np.linspace(-100.0, 100.0, POINTS)

OUTPUT_PREFIX = f"dispersion_run_{RUN}"
# ------------------------------------------


def solve_grid():
    """Select the root with the largest Im(omega) found at each point."""
    omega_grid = np.full(
        (len(KY_VALUES), len(KX_VALUES)),
        complex(np.nan, np.nan),
        dtype=complex,
    )

    for iy, ky in enumerate(KY_VALUES):
        for ix, kx in enumerate(KX_VALUES):
            guesses = []

            # Continuation from previously solved neighbors.
            if ix > 0 and np.isfinite(omega_grid[iy, ix - 1]):
                guesses.append(omega_grid[iy, ix - 1])

            if iy > 0 and np.isfinite(omega_grid[iy - 1, ix]):
                guesses.append(omega_grid[iy - 1, ix])

            roots = find_roots(
                eq_73,
                kx,
                ky,
                RUN,
                guesses=guesses,
                include_default_guesses=True,
            )

            if roots:
                omega_grid[iy, ix] = roots[0]

        print(f"Completed row {iy + 1}/{len(KY_VALUES)}", flush=True)

    return omega_grid


def main():
    omega_grid = solve_grid()

    filename = f"{OUTPUT_PREFIX}_roots.npz"
    np.savez_compressed(
        filename,
        kx=KX_VALUES,
        ky=KY_VALUES,
        omega=omega_grid,
        run=RUN,
        altitude=alt_eff[RUN],
    )

    valid = np.isfinite(omega_grid)

    print(f"\nAccepted grid points: {np.count_nonzero(valid)}")
    print(f"Unstable points: {np.count_nonzero(valid & (omega_grid.imag > 0))}")
    print(f"Damped points: {np.count_nonzero(valid & (omega_grid.imag < 0))}")
    print(f"Missing points: {np.count_nonzero(~valid)}")
    print(f"Saved roots: {filename}")


if __name__ == "__main__":
    main()
