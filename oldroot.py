import os
import matplotlib.pyplot as plt
from scipy.special import wofz
from oldconstants import *

OLD_STUFF_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "old_stuff")
PLOTS_DIR = os.path.join(OLD_STUFF_DIR, "plots")
CSV_DIR = os.path.join(OLD_STUFF_DIR, "csv")
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(CSV_DIR, exist_ok=True)

def fluid_seed_pair_2d(k_x, k_y):
    """Return both fluid-branch roots to seed Newton near both dispersion branches."""
    nu_i_r = nu_i 
    nu_e_r = nu_e 
    kappa_r = kappa 
    k_sq = k_x**2 + k_y**2
    a = 1.0
    b = 1j * nu_i_r + 1j * (m_e / m_i) * nu_e_r * (1 + kappa_r**2)
    c = -C_s_squared * k_sq - 1j * V_d * Omega_i * (kappa_r * k_x - k_y)
    disc = np.sqrt(b**2 - 4 * a * c)
    return ((-b + disc) / (2 * a), (-b - disc) / (2 * a))

def directional_target_root(kx, ky):
    """Select the fluid root corresponding to the unstable wave, respecting the kx sign symmetry."""
    r1, r2 = fluid_seed_pair_2d(kx, ky)
    if kx >= 0:
        return r1 if np.real(r1) > 0 else r2
    else:
        return r1 if np.real(r1) < 0 else r2

def numerical_derivative_2d(func, omega, kx, ky, rel_step=1e-7, min_step=1e-8):
    """Adaptive central difference derivative for 2D complex-valued root finding."""
    scale = max(1.0, abs(omega))
    h = max(min_step, rel_step * scale)
    f_plus = func(omega + h, kx, ky)
    f_minus = func(omega - h, kx, ky)
    return (f_plus - f_minus) / (2 * h)

def analytical_guess(k_x):
    """1D fluid-branch estimate used to seed the 1D Newton-Raphson solver."""
    omega_r = (V_d * k_x) / (1 + Psi_0)
    omega_i = (omega_r**2 - C_s_squared * k_x**2) / (nu_i * (1 + 1 / Psi_0))
    return omega_r + 1j * omega_i

def exact_fluid_root(k_x):
    """1D fluid-branch root (single k_x, fixed k_y) used to seed the 1D Newton-Raphson solver."""
    a = 1.0
    b = 1j * nu_i + 1j * (m_e / m_i) * nu_e * (1 + kappa**2)
    c = -C_s_squared * (k_x**2 + k_y**2) - 1j * V_d * Omega_i * (kappa * k_x - k_y)

    root1 = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
    root2 = (-b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)

    if np.real(root1) > 0:
        return root1
    else:
        return root2

def numerical_derivative(func, omega, k, rel_step=1e-7, min_step=1e-8):
    """Adaptive central difference derivative for 1D complex-valued root finding."""
    scale = max(1.0, abs(omega))
    h = max(min_step, rel_step * scale)
    f_plus = func(omega + h, k)
    f_minus = func(omega - h, k)
    return (f_plus - f_minus) / (2 * h)

def custom_newton_raphson(
    func,
    initial_guess,
    k,
    tol_step=1e-6,
    tol_residual=1e-6,
    max_iter=250,
    max_backtrack=12,
):
    omega_n = initial_guess

    for _ in range(max_iter):
        f_val = func(omega_n, k)
        residual = abs(f_val)

        if residual < tol_residual:
            return omega_n

        df_val = numerical_derivative(func, omega_n, k)

        if abs(df_val) < 1e-14:
            break

        step = -f_val / df_val

        alpha = 1.0
        accepted = False
        for _ in range(max_backtrack):
            omega_trial = omega_n + alpha * step
            if abs(func(omega_trial, k)) < residual:
                accepted = True
                break
            alpha *= 0.5

        if not accepted:
            break

        if abs(alpha * step) < tol_step * (1 + abs(omega_trial)):
            return omega_trial

        omega_n = omega_trial

    return omega_n

def Z(beta):
    # Use asymptotic form for very large |beta| to avoid special-function overflow.
    if abs(beta) > 40:
        return -(1.0 / beta) - (1.0 / (2.0 * beta**3)) - (3.0 / (4.0 * beta**5))

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        z_val = 1j * np.sqrt(np.pi) * wofz(beta)
    return z_val

def Z_tilde(z, growth_rate, nu_in):
    """Piece-wise modification of the Z-function (Eq 72)."""
    # Re(Ri) > 0 is defined as growth_rate > -nu_in
    if growth_rate > -nu_in:
        return Z(z)
    else:
        return -Z(-z)

def Z_prime(beta):
    # Asymptotic derivative of Z(beta): Z'(beta) ~ 1/beta^2 + 3/(2 beta^4) + 15/(4 beta^6)
    if abs(beta) > 40:
        return (1.0 / beta**2) + (3.0 / (2.0 * beta**4)) + (15.0 / (4.0 * beta**6))

    z_beta = Z(beta)

    with np.errstate(over="ignore", invalid="ignore"):
        zp = -2.0 * (1.0 + beta * z_beta)
    return zp

def custom_newton_raphson_2d(
    func,
    initial_guess,
    kx,
    ky,
    tol_step=1e-8,
    tol_residual=1e-6,
    max_iter=250,
    max_backtrack=12,
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

def contour_plot(func, kx_vals, ky_vals, custom_levels, name, new_V_d=None):
    global V_d
    original_V_d = V_d
    func_globals = func.__globals__
    had_func_V_d = 'V_d' in func_globals
    original_func_V_d = func_globals.get('V_d')

    if new_V_d is not None:
        V_d = new_V_d
        if had_func_V_d:
            func_globals['V_d'] = new_V_d

    try:
        return _contour_plot_impl(func, kx_vals, ky_vals, custom_levels, name)
    finally:
        V_d = original_V_d
        if had_func_V_d:
            func_globals['V_d'] = original_func_V_d


def _contour_plot_impl(func, kx_vals, ky_vals, custom_levels, name):
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

            target_root = directional_target_root(kx, ky)

            fluid_r1, fluid_r2 = fluid_seed_pair_2d(kx, ky)
            ref_r1, ref_r2 = fluid_seed_pair_2d(k_mag, 0.0)
            seeds = [fluid_r1, fluid_r2, ref_r1, ref_r2]

            if i > 0:
                seeds.append(root_2d[i - 1, j])
            if j > 0 and root_2d[i, j - 1] != 0j:
                seeds.append(root_2d[i, j - 1])
            if j < len(kx_vals) - 1 and root_2d[i, j + 1] != 0j:
                seeds.append(root_2d[i, j + 1])
            if i > 0 and j > 0 and root_2d[i - 1, j - 1] != 0j:
                seeds.append(root_2d[i - 1, j - 1])

            best_root = target_root
            best_score = np.inf

            for seed in seeds:
                root = custom_newton_raphson_2d(
                func=lambda omega, kx, ky: func(omega, kx, ky),  # The wrapper
                initial_guess=seed,
                kx=kx,
                ky=ky,
                max_iter=NR_MAX_ITER,
                max_backtrack=NR_MAX_BACKTRACK,
                )
                if np.isfinite(root.real) and np.isfinite(root.imag):
                    score = abs(root - target_root) / (1.0 + abs(target_root))
                    if score < best_score:
                        best_score = score
                        best_root = root

            root_2d[i, j] = best_root
            growth_rate_2d[i, j] = np.imag(best_root)

    KX, KY = np.meshgrid(kx_vals, ky_vals)
    growth_rate_data = np.column_stack(
        (KX.ravel(), KY.ravel(), growth_rate_2d.ravel())
    )
    np.savetxt(
        os.path.join(CSV_DIR, f"{name}_growth_rates.csv"),
        growth_rate_data,
        delimiter=",",
        header="kx,ky,growth_rate",
        comments="",
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    cp1 = ax.contourf(KX, KY, growth_rate_2d, levels=custom_levels, cmap='viridis', extend='both') # 'RdBu_r'

    fig.colorbar(cp1, ax=ax)
    ax.set_xlabel('kx (1/m)')
    ax.set_ylabel('ky (1/m)')
    ax.set_title(name)
    ax.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, name), dpi=300)
    plt.show()