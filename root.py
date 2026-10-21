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

    k_sq = k_x**2 + k_y**2
    a = 1.0
    b = 1j * nu_i[run] + 1j * (m_e / m_i) * nu_e[run] * (1 + kappa[run]**2)
    c = -C_s_squared * k_sq - 1j * V_d * Omega_i * (kappa[run] * k_x - k_y)
    disc = np.sqrt(b**2 - 4 * a * c)
    return ((-b + disc) / (2 * a), (-b - disc) / (2 * a))

def numerical_derivative_2d(func, omega, kx, ky, rel_step=1e-7, min_step=1e-8):
    """Adaptive central difference derivative for 2D complex-valued root finding."""
    scale = max(1.0, abs(omega))
    h = max(min_step, rel_step * scale)
    f_plus = func(omega + h, kx, ky)
    f_minus = func(omega - h, kx, ky)
    return (f_plus - f_minus) / (2 * h)

def Z(z):
    return 1j * np.sqrt(np.pi) * wofz(z)

def Z_prime(beta):
    # Asymptotic derivative of Z(beta): Z'(beta) ~ 1/beta^2 + 3/(2 beta^4) + 15/(4 beta^6)
    if abs(beta) > 40:
        return (1.0 / beta**2) + (3.0 / (2.0 * beta**4)) + (15.0 / (4.0 * beta**6))

    z_beta = Z(beta)
    with np.errstate(over="ignore", invalid="ignore"):
        zp = -2.0 * (1.0 + beta * z_beta)
    return zp

def Z_tilde(z, omega_i):
    if np.imag(omega_i) > 0:
        return Z(z)
    elif np.imag(omega_i) < 0:
        return -Z(-z)
    return complex(np.nan, np.nan)

def Z_prime_tilde(z, omega_i):
    if np.imag(omega_i) > 0:
        return Z_prime(z)
    elif np.imag(omega_i) < 0:
        return Z_prime(-z)
    return complex(np.nan, np.nan)

def custom_newton_raphson_2d(
    func, initial_guess, kx, ky, tol_step=1e-8, tol_residual=1e-6, max_iter=250, max_backtrack=12,
    trust_region_scale=2.0,
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
        # Trust-region cap: a near-degenerate Jacobian can otherwise fling the
        # step onto an unrelated, distant root and lose track of the seeded branch.
        max_step = trust_region_scale * max(abs(omega_n), 1.0)
        if abs(step) > max_step:
            step *= max_step / abs(step)
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
    KX, KY = np.meshgrid(kx_vals, ky_vals)
    ny, nx = KX.shape
    growth_rate_2d = np.zeros((ny, nx))
    root_2d = np.zeros((ny, nx), dtype=np.complex128)
    solved = np.zeros((ny, nx), dtype=bool)

    # Process points in order of increasing |k| so continuation only ever
    # extrapolates a short distance from the known-correct k=0 fluid limit,
    # instead of carrying in whatever root was first found at the domain edge.
    order = np.argsort(np.hypot(KX, KY).ravel())
    total = len(order)

    for count, flat_idx in enumerate(order):
        if count % 500 == 0:
            print(f"Point progress: {count}/{total}")

        i, j = np.unravel_index(flat_idx, KX.shape)
        kx, ky = KX[i, j], KY[i, j]

        fluid_r1, fluid_r2 = fluid_seed_pair_2d(kx, ky, run)
        seeds = [fluid_r1, fluid_r2]

        # Seed continuation only from already-solved 8-connected neighbors
        # (no more "reference at (k_mag,0)" seeds -- those belong to a
        # different (kx,ky) location's fluid solution and can converge to a
        # spurious, unrelated root even at very small |k|).
        neighbors = []
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == 0 and dj == 0:
                    continue
                ni, nj = i + di, j + dj
                if 0 <= ni < ny and 0 <= nj < nx and solved[ni, nj]:
                    dist = np.hypot(KX[ni, nj] - kx, KY[ni, nj] - ky)
                    neighbors.append((dist, root_2d[ni, nj]))

        seeds.extend(r for _, r in neighbors)

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
            if neighbors:
                # Anchor on the nearest (by physical kx,ky distance) already-solved
                # neighbor to stay on the same branch instead of jumping to a
                # faster-growing but unrelated root.
                anchor_root = min(neighbors, key=lambda t: t[0])[1]
                best_root = min(valid_roots, key=lambda r: abs(r - anchor_root))
            else:
                best_root = max(valid_roots, key=lambda r: np.imag(r))
        else:
            best_root = 0j

        root_2d[i, j] = best_root
        growth_rate_2d[i, j] = np.imag(best_root)
        solved[i, j] = True

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