import os
import matplotlib.pyplot as plt
from scipy.special import wofz
import numpy as np
from constants import *

PLOTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")
CSV_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "csv")
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(CSV_DIR, exist_ok=True)

def fluid_seed_pair_2d(k_x, k_y, run=None):
    """Return both fluid-branch roots to seed Newton near both dispersion branches."""
    nu_i_r = nu_i if run is None else nu_i[run]
    nu_e_r = nu_e if run is None else nu_e[run]
    kappa_r = kappa if run is None else kappa_e[run]
    k_sq = k_x**2 + k_y**2
    a = 1.0
    b = 1j * nu_i_r + 1j * (m_e / m_i) * nu_e_r * (1 + kappa_r**2)
    c = -C_s_squared * k_sq - 1j * V_d * Omega_i * (kappa_r * k_x - k_y)
    disc = np.sqrt(b**2 - 4 * a * c)
    return ((-b + disc) / (2 * a), (-b - disc) / (2 * a))

def numerical_derivative_2d(func, omega, kx, ky, rel_step=1e-7, min_step=1e-8):
    """Adaptive central difference derivative for 2D complex-valued root finding."""
    scale = max(1.0, abs(omega))
    h = max(min_step, rel_step * scale)
    f_plus = func(omega + h, kx, ky)
    f_minus = func(omega - h, kx, ky)
    return (f_plus - f_minus) / (2 * h)

def Z(beta):
    if abs(beta) > 40:
        return -(1.0 / beta) - (1.0 / (2.0 * beta**3)) - (3.0 / (4.0 * beta**5))
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        z_val = 1j * np.sqrt(np.pi) * wofz(beta)
    return z_val

def Z_tilde(z, growth_rate, nu_in):
    if growth_rate > -nu_in:
        return Z(z)
    else:
        return -Z(-z)

def custom_newton_raphson_2d(
    func, initial_guess, kx, ky, tol_step=1e-8, tol_residual=1e-6, max_iter=250, max_backtrack=12
):
    omega_n = initial_guess
    for _ in range(max_iter):
        f_val = func(omega_n, kx, ky)
        residual = abs(f_val)
        if residual < tol_residual:
            return omega_n
        df_val = numerical_derivative_2d(func, omega_n, kx, ky)
        if abs(df_val) < 1e-14:
            break
        step = -f_val / df_val
        alpha = 1.0
        accepted = False
        for _ in range(max_backtrack):
            omega_trial = omega_n + alpha * step
            f_trial = func(omega_trial, kx, ky)
            if abs(f_trial) < residual:
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            break
        if abs(alpha * step) < tol_step * (1 + abs(omega_trial)):
            return omega_trial
        omega_n = omega_trial
    return omega_n

def contour_plot(func, kx_vals, ky_vals, name, pts, run=None, custom_levels=None, positive_only=False):
    growth_rate_2d = np.zeros((len(ky_vals), len(kx_vals)))
    root_2d = np.zeros((len(ky_vals), len(kx_vals)), dtype=np.complex128)

    center_j = np.argmin(np.abs(kx_vals))

    for i, ky in enumerate(ky_vals):
        if i % 10 == 0:
            print(f"Row progress: {i}/{len(ky_vals)}")

        j_order = list(range(center_j, len(kx_vals))) + list(range(center_j - 1, -1, -1))

        for j in j_order:
            kx = kx_vals[j]
            k_mag = np.hypot(kx, ky)

            fluid_r1, fluid_r2 = fluid_seed_pair_2d(kx, ky, run)
            ref_r1, ref_r2 = fluid_seed_pair_2d(k_mag, 0.0, run)
            seeds = [fluid_r1, fluid_r2, ref_r1, ref_r2]

            if i > 0:
                seeds.append(root_2d[i - 1, j])
            if j > 0 and root_2d[i, j - 1] != 0j:
                seeds.append(root_2d[i, j - 1])
            if j < len(kx_vals) - 1 and root_2d[i, j + 1] != 0j:
                seeds.append(root_2d[i, j + 1])
            if i > 0 and j > 0 and root_2d[i - 1, j - 1] != 0j:
                seeds.append(root_2d[i - 1, j - 1])

            valid_roots = []

            for seed in seeds:
                root = custom_newton_raphson_2d(
                    func=lambda omega, kx, ky: func(omega, kx, ky),
                    initial_guess=seed,
                    kx=kx,
                    ky=ky,
                    max_iter=NR_MAX_ITER,
                    max_backtrack=NR_MAX_BACKTRACK,
                )

                if np.isfinite(root.real) and np.isfinite(root.imag):
                    if abs(func(root, kx, ky)) < 1e-3:
                        valid_roots.append(root)

            if valid_roots:
                best_root = max(valid_roots, key=lambda r: np.imag(r))
            else:
                best_root = 0j

            root_2d[i, j] = best_root
            growth_rate_2d[i, j] = np.imag(best_root)

    KX, KY = np.meshgrid(kx_vals, ky_vals)
    growth_rate_data = np.column_stack((KX.ravel(), KY.ravel(), growth_rate_2d.ravel()))
    #np.savetxt(os.path.join(CSV_DIR, f"{name}_growth_rates.csv"), growth_rate_data, delimiter=",", header="kx,ky,growth_rate", comments="")
    #np.savetxt(os.path.join(CSV_DIR, f"{name}_a0_zero_growth_rates.csv"), growth_rate_data, delimiter=",", header="kx,ky,growth_rate", comments="")

    fig, ax = plt.subplots(figsize=(8, 6))

    if positive_only:
        growth_rate_2d = np.where(growth_rate_2d <= 0, np.nan, growth_rate_2d)
        if custom_levels is not None:
            levels = custom_levels
        else:
            hi = np.nanmax(growth_rate_2d) if np.nanmax(growth_rate_2d) > 0 else 8000
            levels = np.linspace(0, hi, 60)
        cp1 = ax.contourf(KX, KY, growth_rate_2d, levels=levels, cmap='viridis', extend='max')
    else:
        # Same rendering convention as eq 2.46: show negative growth too, fixed-or-auto levels, extend both ends.
        if custom_levels is not None:
            levels = custom_levels
        else:
            lo = min(0.0, growth_rate_2d.min())
            hi = growth_rate_2d.max() if growth_rate_2d.max() > 0 else 8000
            levels = np.linspace(lo, hi, 60)
        cp1 = ax.contourf(KX, KY, growth_rate_2d, levels=levels, cmap='viridis', extend='both')

    fig.colorbar(cp1, ax=ax)
    ax.set_xlabel('kx (1/m)')
    ax.set_ylabel('ky (1/m)')
    ax.set_title(fr'Effective altitude {name} km')
    ax.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, f"TEST {name} km"), dpi=300)
    plt.show()