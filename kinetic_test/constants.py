"""
TEST constants for the kinetic-vs-fluid comparison.

Lives in its own folder (kinetic_test/) next to a COPY of root.py, so that
`from constants import *` inside root.py picks up THIS file and your real
constants.py is never touched.

!!! Every number below marked (PLACEHOLDER) is a typical E-region value I made up.
!!! Replace them with the values from your real constants.py before drawing any
!!! conclusions.  Only the structure matters for debugging the equation.

Temperatures are in JOULES (k_B * T[K]) because the solver uses sqrt(T/m).
All per-run quantities are length-1 lists, so use run = 0 everywhere.
"""
import numpy as np

# ---------------------------------------------------------------- physical constants
e      = 1.602176634e-19
m_e    = 9.1093837015e-31

k_B    = 1.380649e-23
eps0   = 8.8541878128e-12

# ---------------------------------------------------------------- background
m_i    = 5e-26          # NO+/O2+
m_n    = m_i/1250          # neutrals
T_n    = 300 * k_B         # neutral temperature [J]
n0     = 4e8              # plasma density [m^-3]
B_0    = 50000e-9              # [T]
V_d = 2100
E_0    =  V_d * B_0                # [V/m]

# collision frequencies, one entry per run (PLACEHOLDER)
nu_i = np.array([488])         # Ion-neutral collision freq (s^-1)
nu_e = np.array([4148])           # electron-neutral [1/s]

Omega_i  = e * B_0 / m_i
Omega_e  = e * B_0 / m_e
kappa_i  = [Omega_i / nu_i[0]]
kappa_e  = [Omega_e / nu_e[0]]
kappa_e = np.array([48])            # frequency ratio
kappa_i = np.array([0.33])

# heating "reduced masses" (CHECK against your own definition)
m_n_e = m_e * m_n / (m_e + m_n)
m_n_i = m_i * m_n / (m_i + m_n)

anisotropy_ratio = 0.5      # T_x / T_y of the ion bi-Maxwellian
gamma = 1.0                  # electron polytropic index

# ---------------------------------------------------------------- derived (cold reference, T_e = T_n)
V_Tn            = np.sqrt(T_n / m_i)
beta_T          = 1.0                                   # T_n / T_e with T_e = T_n
lambda_D_sq     = eps0 * T_n / (e**2 * n0)
omega_pe_squared = e**2 * n0 / (eps0 * m_e)

# names that root.py's fluid_seed_pair_2d / contour_plot expect
C_s_squared = (gamma * T_n + T_n) / m_i
V_d         = E_0 / B_0
kappa       = kappa_e[0]
NR_MAX_ITER     = 250
NR_MAX_BACKTRACK = 12

# ---------------------------------------------------------------- TEST SWITCHES
# Edit these one at a time.  1 = physics on, 0 = off.
FLAGS = {
    "HEAT_ELECTRONS": 0,          # T_e = T_n + (m_n_e/3) V_de^2 ;  0 -> T_e = T_n
    "HEAT_IONS": 0,               # T_eff_i from frictional heating ; 0 -> T_n
    "ION_ANISOTROPY": 0,          # use anisotropy_ratio ; 0 -> T_x = T_y
    "ION_FIELD_TERMS": 0,         # the i k.a_i0 terms in den1, den2
    "ION_DRIFT": 0,               # U0 in the second Z-function
    "ELECTRON_SIGN": +1,          # +1 = derived sign, -1 = your original code
    "CONSISTENT_NORMALIZATION": 1 # recompute lambda_D^2 and beta_T from the local T_e
}
