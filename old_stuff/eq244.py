import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from oldroot import *

def eq_244(omega, k_x):
    beta = (omega + 1j * nu_i) / (np.sqrt(2) * V_th * k_x)

    term1 = omega_pe_squared / (-1j * omega * nu_e * (1 + kappa**2) + (gamma * T_e / m_e) * k_x**2 + 1j * nu_e * kappa**2 * V_d * k_x)
    term2 = (omega_pi_squared * Z_prime(beta)) / (2 * V_th**2 * k_x**2 + np.sqrt(2) * 1j * nu_i * V_th * k_x * Z(beta))

    return term1 - term2 - 1.0

k_list_growth = np.logspace(-1, 2, 500)
k_list_freq = np.linspace(0.1, 30, 200)

omega_list_kinetic_growth = []
prev_root = None

for k in k_list_growth:
    if prev_root is None:
        start_guess = exact_fluid_root(k)
    else:
        start_guess = prev_root

    root = custom_newton_raphson(
        func=eq_244,
        initial_guess=start_guess,
        k=k,
    )
    omega_list_kinetic_growth.append(root)
    prev_root = root

omega_list_kinetic_freq = []
prev_root = None

for k in k_list_freq:
    if prev_root is None:
        start_guess = exact_fluid_root(k)
    else:
        start_guess = prev_root

    root = custom_newton_raphson(
            func=eq_244,
            initial_guess=start_guess,
            k=k,
        )
    omega_list_kinetic_freq.append(root)
    prev_root = root

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(k_list_growth, np.imag(omega_list_kinetic_growth), label="Kinetic")
ax1.set_xscale("log")
ax1.set_xlabel("k (1/m)")
ax1.set_ylabel("growth rate (1/s)")
ax1.set_xlim(0.1, 100)
ax1.set_ylim(-100, 400)
ax1.grid(True)
ax1.legend()

ax2.plot(k_list_freq, np.real(omega_list_kinetic_freq), label="Kinetic")
ax2.set_xlabel("k (1/m)")
ax2.set_ylabel("freq (rad/s)")
ax2.set_xlim(0, 30)
ax2.set_ylim(0, 12000)
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(__file__), "plots", "sol_eq244.png"))
plt.show()