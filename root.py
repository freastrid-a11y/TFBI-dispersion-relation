import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import root
from scipy.special import wofz

# Exposes the constants to eq73.py through its existing:
#     from root import *
from constants import *


def Z_tilde(z, omega_i=None):
    """
    ASSUMPTION: Z_tilde is the standard plasma dispersion function.

        Z(z) = i * sqrt(pi) * wofz(z)

    omega_i is accepted for compatibility with eq73.py, but is not used.
    Replace this function if your model defines a modified Z_tilde.
    """
    with np.errstate(over="ignore", invalid="ignore"):
        return 1j * np.sqrt(np.pi) * wofz(z)


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

    # eq_73 already multiplies its Poisson sum by k^2 * lambda_D_sq.
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


def contour_plot(
    dispersion,
    kx_vals,
    ky_vals,
    altitude,
    pts=None,
    run=0,
    positive_only=True,
):
    """
    Plot the largest growth rate found at each wave-vector grid point.

    Assumes perturbations vary as exp(i*k·r - i*omega*t), so Im(omega) > 0
    indicates growth.

    pts is accepted for compatibility with your existing eq73.py.
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

        print(
            f"Altitude {altitude:g} km: "
            f"row {iy + 1}/{len(ky_vals)}"
        )

    growth = omega_grid.imag.copy()

    if positive_only:
        growth[growth <= 0.0] = np.nan

    growth = np.ma.masked_invalid(growth)

    fig, ax = plt.subplots(figsize=(8, 6))

    if growth.count() == 0:
        ax.text(
            0.5,
            0.5,
            "No accepted unstable roots found"
            if positive_only else "No accepted roots found",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    else:
        mesh = ax.pcolormesh(
            kx_vals,
            ky_vals,
            growth,
            shading="auto",
            cmap="viridis",
        )

        fig.colorbar(
            mesh,
            ax=ax,
            label=r"Growth rate $\mathrm{Im}(\omega)$ [s$^{-1}$]",
        )

    ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
    ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    ax.set_title(f"Eq. 73 — altitude {altitude:g} km")
    plt.savefig("104.png")

    fig.tight_layout()
    plt.show()

    return omega_grid