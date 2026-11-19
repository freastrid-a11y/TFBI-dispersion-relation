import numpy as np
from scipy.optimize import root
from scipy.special import wofz
from constants import *


def Z_tilde(z, omega_i):
    """
    Modified plasma dispersion function from manuscript Eq. (74).

    The same omega_i determines the branch for both ion arguments.
    The manuscript does not specify the value at Im(omega_i) = 0.
    """
    with np.errstate(over="ignore", invalid="ignore"):
        if np.imag(omega_i) > 0.0:
            # Z(z)
            return 1j * np.sqrt(np.pi) * wofz(z)

        if np.imag(omega_i) < 0.0:
            # -Z(-z)
            return -1j * np.sqrt(np.pi) * wofz(-z)

        return complex(np.nan, np.nan)


def frequency_scale(k_x, k_y, run):
    """Characteristic frequency used to nondimensionalize omega."""
    k_mag = np.hypot(k_x, k_y)

    return float(max(
        1.0,
        nu_i[run],
        k_mag * abs(V_d),
        k_mag * V_Tx[run],
        k_mag * V_Ty[run],
    ))


def default_guesses(k_x, k_y, run, scale):
    """Initial guesses spanning several phase velocities and growth rates."""
    V_de_x = V_d * kappa_e[run]**2 / (1.0 + kappa_e[run]**2)
    V_de_y = -V_d * kappa_e[run] / (1.0 + kappa_e[run]**2)

    omega_e_drift = k_x * V_de_x + k_y * V_de_y
    omega_i_drift = k_y * U0_y[run]

    # Approximate low-frequency drift-wave phase frequency.
    omega_drift = (
        omega_e_drift + Psi_0[run] * omega_i_drift
    ) / (1.0 + Psi_0[run])

    real_guesses = np.unique([
        0.0,
        omega_drift,
        omega_e_drift,
        omega_i_drift,
        -0.5 * scale,
        0.5 * scale,
    ])

    imag_guesses = [
        0.05 * scale,
        -0.05 * scale,
        -0.5 * scale,
    ]

    return [
        complex(wr, wi)
        for wr in real_guesses
        for wi in imag_guesses
    ]


def find_roots(
    dispersion,
    k_x,
    k_y,
    run,
    guesses=None,
    include_default_guesses=True,
    residual_tol=1e-8,
    correction_tol=1e-6,
    max_frequency_factor=30.0,
):
    """
    Find candidate roots of dispersion(omega, k_x, k_y, run) = 0.

    Returns
    -------
    list of complex
        Distinct accepted roots, sorted by decreasing imaginary part.

    Notes
    -----
    This is a local, multistart solver, not a proof that all roots were found.
    max_frequency_factor restricts accepted roots to the searched frequency
    range and helps reject spurious convergence at extremely large |omega|.
    """
    if k_x == 0.0 and k_y == 0.0:
        # The equation regularizes k^2, but den2 still vanishes at k = 0.
        return []

    scale = frequency_scale(k_x, k_y, run)

    seeds = [] if guesses is None else list(guesses)

    if include_default_guesses:
        seeds.extend(default_guesses(k_x, k_y, run, scale))

    if not seeds:
        seeds = default_guesses(k_x, k_y, run, scale)

    # dispersion_rel already multiplies its Poisson sum by k^2 * lambda_D_sq.
    # Use an additional FIXED normalization so the solver tolerance remains
    # meaningful across different k values.
    k_sq = max(k_x**2 + k_y**2, 1e-12)

    electron_reference = omega_pe_squared / (
        nu_e[run] * (1.0 + kappa_e[run]**2) * scale
    )

    residual_scale = (
        k_sq * lambda_D_sq * max(1.0, electron_reference)
    )

    def evaluate(omega):
        with np.errstate(
            over="ignore",
            divide="ignore",
            invalid="ignore",
            under="ignore",
        ):
            value = dispersion(omega, k_x, k_y, run)

        return complex(value) / residual_scale

    def residual(x):
        omega = scale * complex(x[0], x[1])
        value = evaluate(omega)

        if not np.isfinite(value):
            # Keep invalid special-function evaluations out of the solver.
            return np.array([1e20, 1e20])

        return np.array([value.real, value.imag])

    accepted = []

    for guess in seeds:
        if not np.isfinite(guess):
            continue

        x0 = np.array([guess.real, guess.imag]) / scale

        try:
            solution = root(
                residual,
                x0,
                method="hybr",
                options={
                    "xtol": 1e-10,
                    "maxfev": 600,
                },
            )
        except (ValueError, FloatingPointError, OverflowError):
            continue

        omega = scale * complex(solution.x[0], solution.x[1])

        if not np.isfinite(omega):
            continue

        if abs(omega) > max_frequency_factor * scale:
            continue

        value = evaluate(omega)

        # Independently check the residual; solver status alone is insufficient.
        if not np.isfinite(value) or abs(value) > residual_tol:
            continue

        # Check the estimated remaining Newton correction as well.
        # This helps reject a small residual caused by a nearly flat function.
        h = 1e-5 * scale
        derivative = (
            evaluate(omega + h) - evaluate(omega - h)
        ) / (2.0 * h)

        if not np.isfinite(derivative) or abs(derivative) == 0.0:
            continue

        correction = abs(value / derivative)

        if correction > correction_tol * scale:
            continue

        # Merge roots reached from different initial guesses.
        duplicate = any(
            abs(omega - previous)
            <= 1e-5 * max(scale, abs(omega), abs(previous))
            for previous in accepted
        )

        if not duplicate:
            accepted.append(omega)

    return sorted(accepted, key=lambda w: w.imag, reverse=True)


def _follow_root(dispersion, run, omega, k_from, k_to, substeps=8):
    """
    Follow a root from k_from to k_to in small k-steps, re-solving with
    find_roots at each step. Returns the root at k_to, or None if it is lost.
    """
    for t in np.linspace(0.0, 1.0, substeps + 1)[1:]:
        kx = (1.0 - t) * k_from[0] + t * k_to[0]
        ky = (1.0 - t) * k_from[1] + t * k_to[1]

        roots = find_roots(
            dispersion, kx, ky, run,
            guesses=[omega], include_default_guesses=False,
        )
        if not roots:
            return None

        # Stay on the same branch: take the root closest to the previous one.
        nearest = min(roots, key=lambda w: abs(w - omega))
        if abs(nearest - omega) > 0.25 * max(abs(omega), 1.0):
            return None
        omega = nearest

    return omega


def solve_grid(dispersion, kx_vals, ky_vals, run=0, max_passes=20, substeps=8):
    """
    Solve for the accepted root with the largest Im(omega) at each
    wave-vector grid point, continuing from already-solved neighbors.

    Assumes perturbations vary as exp(i*k·r - i*omega*t), so Im(omega) > 0
    indicates growth. Returns the complex omega grid with shape
    (len(ky_vals), len(kx_vals)).

    After the row-by-row sweep, the multistart search can still have returned
    a damped root where a growing one exists (e.g. the rectangular block-out
    seen in 104 km double-E runs). Neighbor passes then re-seed every stable
    point that touches an unstable one:

      * from each unstable neighbor (all 8), by following that neighbor's
        root to this point in `substeps` small k-steps;
      * from the mirror point, omega(-k) = -conj(omega(k)), when the grid
        is symmetric about k = 0.

    A point is updated only if this gives an accepted root (same find_roots
    tests) with a larger Im(omega). Passes repeat until nothing changes, at
    most max_passes times.
    """
    kx_vals = np.asarray(kx_vals, dtype=float)
    ky_vals = np.asarray(ky_vals, dtype=float)

    ny, nx = len(ky_vals), len(kx_vals)
    omega_grid = np.full((ny, nx), np.nan + 1j * np.nan, dtype=complex)

    # --- Sweep: row by row, seeded from the left and lower neighbors ---
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

        print(f"run {run}: row {iy + 1}/{ny}", flush=True)

    # --- Neighbor passes: re-seed stable points from unstable neighbors ---
    symmetric = (np.allclose(kx_vals, -kx_vals[::-1])
                 and np.allclose(ky_vals, -ky_vals[::-1]))

    for n_pass in range(max_passes):
        unstable = np.isfinite(omega_grid) & (omega_grid.imag > 0.0)
        snapshot = omega_grid.copy()
        improved = 0

        for iy in range(ny):
            for ix in range(nx):
                if unstable[iy, ix]:
                    continue

                neighbors = [
                    (iy + dy, ix + dx)
                    for dy in (-1, 0, 1) for dx in (-1, 0, 1)
                    if (dy or dx)
                    and 0 <= iy + dy < ny and 0 <= ix + dx < nx
                    and unstable[iy + dy, ix + dx]
                ]
                mirror = symmetric and unstable[ny - 1 - iy, nx - 1 - ix]

                if not neighbors and not mirror:
                    continue

                kx, ky = kx_vals[ix], ky_vals[iy]
                candidates = []

                if mirror:
                    seed = -np.conj(snapshot[ny - 1 - iy, nx - 1 - ix])
                    candidates += find_roots(
                        dispersion, kx, ky, run,
                        guesses=[seed], include_default_guesses=False,
                    )

                for jy, jx in neighbors:
                    followed = _follow_root(
                        dispersion, run, snapshot[jy, jx],
                        (kx_vals[jx], ky_vals[jy]), (kx, ky), substeps,
                    )
                    if followed is not None:
                        candidates += find_roots(
                            dispersion, kx, ky, run,
                            guesses=[followed], include_default_guesses=False,
                        )

                if not candidates:
                    continue

                best = max(candidates, key=lambda w: w.imag)
                current = omega_grid[iy, ix].imag if np.isfinite(omega_grid[iy, ix]) else -np.inf

                if best.imag > current + 1e-9 * abs(best):
                    omega_grid[iy, ix] = best
                    improved += 1

        print(f"run {run}: neighbor pass {n_pass + 1}: {improved} points improved", flush=True)

        if improved == 0:
            break

    return omega_grid