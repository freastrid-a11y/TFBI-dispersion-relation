import numpy as np
from scipy.optimize import root
from scipy.special import wofz
from constants import *


def Z_tilde(z, omega_i):
    """
    Modified plasma dispersion function from Eq. (74).
    """
    with np.errstate(over="ignore", invalid="ignore"):
        if np.imag(omega_i) > 0.0:
            # Z(z)
            return 1j * np.sqrt(np.pi) * wofz(z)

        if np.imag(omega_i) < 0.0:
            # -Z(-z)
            return -1j * np.sqrt(np.pi) * wofz(-z)

        return complex(np.nan, np.nan)


def frequency_scale(k_x, k_y, run): # non-dimensionalizing omega
    k_mag = np.hypot(k_x, k_y)

    return float(max(
        1.0,
        nu_i[run],
        k_mag * abs(V_d),
        k_mag * V_Tx[run],
        k_mag * V_Ty[run],
    ))


def default_guesses(k_x, k_y, run, scale): # generate guesses
    # Electron drift
    V_de_x = V_d * kappa_e[run]**2 / (1.0 + kappa_e[run]**2)
    V_de_y = -V_d * kappa_e[run] / (1.0 + kappa_e[run]**2)

    # Doppler frequencies
    omega_e_drift = k_x * V_de_x + k_y * V_de_y
    omega_i_drift = k_y * U0_y[run]

    omega_drift = ( # real frequency of the linear fluid Farley-Bunemann instability
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
    residual_tol=1e-8,
    correction_tol=1e-6,
    max_frequency_factor=30.0,
):
    """
    Find the roots omega of dispersion(omega, k_x, k_y, run) = 0.

    Tries the given guesses plus default_guesses and keeps every
    distinct root that passes the checks. Returns a list sorted with the
    most unstable root (largest Im omega) first.

    Only finds roots near the guesses, so it can miss some. Roots with
    |omega| > max_frequency_factor * scale are thrown away.
    """
    if k_x == 0.0 and k_y == 0.0:
        return []  # the ion term divides by zero at k = 0

    scale = frequency_scale(k_x, k_y, run)

    seeds = [] if guesses is None else list(guesses)
    seeds.extend(default_guesses(k_x, k_y, run, scale))

    # Divide D by an estimate of its typical size, so that one residual_tol
    # can be used at every k. The estimate is built for the fluid-electron
    # form, D = k^2 lambda_D^2 (1 + chi_e + chi_i): electron_reference is
    # roughly |chi_e|, and k^2 lambda_D^2 multiplies everything.
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
            # Return a huge value so the solver moves away from this omega.
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

        # Check |D| ourselves; the solver can stop without really reaching 0.
        if not np.isfinite(value) or abs(value) > residual_tol:
            continue

        # Estimate the next Newton step, |D / D'|. If it is still large, the
        # small |D| came from a flat region of D, not from a real root.
        h = 1e-5 * scale
        derivative = (
            evaluate(omega + h) - evaluate(omega - h)
        ) / (2.0 * h)

        if not np.isfinite(derivative) or abs(derivative) == 0.0:
            continue

        correction = abs(value / derivative)

        if correction > correction_tol * scale:
            continue

        # Different guesses often land on the same root; keep it only once.
        duplicate = any(
            abs(omega - previous)
            <= 1e-5 * max(scale, abs(omega), abs(previous))
            for previous in accepted
        )

        if not duplicate:
            accepted.append(omega)

    return sorted(accepted, key=lambda w: w.imag, reverse=True)


def solve_grid(dispersion, kx_vals, ky_vals, run=0):
    """
    Find the most unstable root (largest Im omega) at every (k_x, k_y).

    Waves vary as exp(i k.r - i omega t), so Im(omega) > 0 means growth.
    Each point also uses the roots already found at its left and lower
    neighbours as starting guesses, which keeps the solution on the same
    branch from point to point.

    Returns a complex array of shape (len(ky_vals), len(kx_vals)),
    NaN where no root was found.
    """
    kx_vals = np.asarray(kx_vals, dtype=float)
    ky_vals = np.asarray(ky_vals, dtype=float)

    shape = (len(ky_vals), len(kx_vals))
    omega_grid = np.full(shape, np.nan + 1j * np.nan, dtype=complex)

    for iy, ky in enumerate(ky_vals):
        for ix, kx in enumerate(kx_vals):
            continuation_seeds = []

            # Use the neighbours' roots (left and below) as extra guesses.
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
            )

            if roots:
                omega_grid[iy, ix] = roots[0]

        print(f"run {run}: row {iy + 1}/{len(ky_vals)}", flush=True)

    return omega_grid