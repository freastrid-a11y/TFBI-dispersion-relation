import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from oldroot import *

def eq_223(omega, k_x):
    term1 = omega**2 + 1j * omega * nu_i - C_s_squared * (k_x**2 + k_y**2)
    term2 = -1j * V_d * Omega_i * (kappa * k_x - k_y)
    term3 = (m_e / m_i) * nu_e * (1 + kappa**2) * 1j * omega
    return term1 + term2 + term3

k_list_growth = np.logspace(-1, 2, 500)
k_list_freq = np.linspace(0.1, 30, 200)

omega_list_growth = []
prev_root_growth = None
for k in k_list_growth:
    guess = analytical_guess(k)
    start_guess = prev_root_growth if prev_root_growth is not None else guess
    root = custom_newton_raphson(eq_223, start_guess, k)
    omega_list_growth.append(root)
    prev_root_growth = root

omega_list_freq = []
prev_root_freq = None
for k in k_list_freq:
    guess = analytical_guess(k)
    start_guess = prev_root_freq if prev_root_freq is not None else guess
    root = custom_newton_raphson(eq_223, start_guess, k)
    omega_list_freq.append(root)
    prev_root_freq = root

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

growth_vals = np.imag(omega_list_growth)
freq_vals = np.real(omega_list_freq)

ax1.plot(k_list_growth, growth_vals, linestyle='--', label='Fluid')
ax1.set_xscale('log')
ax1.set_xlim(0.1, 100)
ax1.set_ylim(-100, 400)
ax1.set_xlabel('k (1/m)')
ax1.set_ylabel('growth rate (1/s)')
ax1.grid(True)
ax1.legend()

ax2.plot(k_list_freq, freq_vals, linestyle='--', label='Fluid')
ax2.set_xlim(0, 30)
ax2.set_ylim(0, 12000)
ax2.set_xlabel('k (1/m)')
ax2.set_ylabel('freq (rad/s)')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(__file__), "plots", "sol_eq223.png"))
plt.show()