# `constants.py` — Physical Constants & Altitude Profiles

## Overview

This module defines every physical constant and derived parameter used by the solver. It is split into three groups:

1. **Fundamental constants and plasma parameters** — masses, charge, permittivity, magnetic field, drift speed, temperatures, densities, plasma frequencies, Debye length.
2. **Altitude-dependent arrays** — an analytic fit $\psi(h)$ to the E-region collision profile generates, for each effective altitude in `alt_eff`, the magnetization parameters $\kappa_{e,i}$, collision frequencies $\nu_{e,i}$, and the driving parameter $\Psi_0$. **Every "run" index used elsewhere in the code (`--run N`) is an index into these arrays.**
3. **Drift/heating-derived quantities** — the ion Pedersen drift $\mathbf{U}_0$, ion acceleration $\mathbf{a}_{i0}$, frictionally heated effective temperature $T_\mathrm{eff}$, and anisotropic thermal velocities $V_{Tx}, V_{Ty}$.

!!! warning "Temperatures are in Joules"
    All temperatures (`T_n`, `T_e`, `T_eff`, `T_x`, `T_y`) are stored as energies, i.e. already multiplied by $k_B$. Thermal velocities are then simply $\sqrt{T/m}$.

### Key physics

- **Electron heating.** The electron temperature includes ohmic heating balanced by inelastic losses:
  $T_e = T_n + \dfrac{2 m_e E_0^2}{3\,\delta_{en} B_0^2}$, with energy-loss fraction $\delta_{en} = 3\times10^{-3}$. (The module prints $T_e/T_n$ at import time as a sanity check.)
- **Collision profile fit.** $\psi(h)$ is a two-exponential analytic fit smoothly switching around 120 km; from it,
  $\kappa_i = \Theta_0/\sqrt{\psi}$, $\kappa_e = 1/(\Theta_0\sqrt{\psi})$, $\nu_{i,e} = \Omega_{i,e}/\kappa_{i,e}$, and $\Psi_0 = 1/(\kappa_i\kappa_e)$.
- **Frictional ion heating.** $T_\mathrm{eff} = T_n + m_i U_0^2/3$, partitioned by `anisotropy_ratio` into $T_x = r\,T_y$ and $T_y = 3T_\mathrm{eff}/(1+2r)$ (currently $r = 1$: isotropic).
- **Geometry.** $V_d$ (electron drift) is along $x$; $E_0$, $\mathbf{U}_0$, and $\mathbf{a}_{i0}$ are along $y$.

---

## API Reference (module-level constants)

### Scalars

| Name | Value / formula | Units | Description |
|---|---|---|---|
| `m_i`, `m_e` | `5e-26`, `9.1094e-31` | kg | Ion and electron masses |
| `m_n_e`, `m_n_i` | `1e-27`, `5e-26` | kg | Neutral masses for e/i collisions |
| `e` | `1.6e-19` | C | Elementary charge |
| `eps_0` | `8.854e-12` | F/m | Vacuum permittivity |
| `k_B` | `1.380649e-23` | J/K | Boltzmann constant |
| `delta_en` | `3e-3` | — | Electron energy-loss fraction |
| `B_0` | `5e-5` | T | Background magnetic field |
| `V_d` | `2100` | m/s | Electron $\mathbf{E}\times\mathbf{B}$ drift speed |
| `E_0` | `V_d * B_0` | V/m | Background electric field |
| `gamma` | `1.0` | — | Electron polytropic index (isothermal) |
| `T_n` | `300 k_B` | J | Neutral/ion background temperature |
| `T_e` | see above | J | Heated electron temperature |
| `beta_T` | `T_n / T_e` | — | Temperature ratio |
| `C_s_squared` | `(T_n + γT_e)/m_i` | m²/s² | Ion-acoustic speed squared |
| `V_th`, `V_Tn` | `sqrt(T_n/m_i)` | m/s | Ion/neutral thermal speed |
| `Omega_i`, `Omega_e` | `eB_0/m` | rad/s | Gyrofrequencies |
| `n_0` | `1e11` | m⁻³ | Plasma density |
| `omega_pe_squared`, `omega_pi_squared` | `n_0 e²/(ε_0 m)` | rad²/s² | Plasma frequencies squared |
| `lambda_D_sq` | `ε_0 T_e/(n_0 e²)` | m² | Debye length squared |
| `Theta_0` | `1.35e-2` | — | Mass-collision ratio constant |
| `anisotropy_ratio` | `1` | — | Ion temperature anisotropy $T_x/T_y$ |
| `a_i0_mag`, `a_i0_x`, `a_i0_y` | `eE_0/m_i`; (0, mag) | m/s² | Ion acceleration vector |
| `NR_MAX_ITER`, `NR_MAX_BACKTRACK` | `140`, `10` | — | Legacy Newton-solver limits |

### Altitude-indexed arrays (length 8, indexed by `run`)

| Name | Formula | Units | Description |
|---|---|---|---|
| `alt_eff` | `[135, 129, 122, 116, 110, 104, 97, 91]` | km | Effective altitudes |
| `psi_h` | analytic fit | — | Collision profile $\psi(h)$ |
| `kappa_i`, `kappa_e` (`kappa`) | `Θ_0/√ψ`, `1/(Θ_0√ψ)` | — | Magnetization parameters |
| `nu_i`, `nu_e` | `Ω/κ` | s⁻¹ | Collision frequencies |
| `Psi_0` | `1/(κ_i κ_e)` | — | FBI driving parameter |
| `U0_mag`, `U0_x`, `U0_y` | `eE_0/(m_i ν_i)`; (0, mag) | m/s | Ion Pedersen drift |
| `T_eff`, `T_x`, `T_y`, `T_eff_i` | see Overview | J | Heated / anisotropic ion temperatures |
| `V_Tx`, `V_Ty` | `sqrt(T/m_i)` | m/s | Anisotropic ion thermal speeds |

### Run-index ↔ altitude map

| `run` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| altitude (km) | 135 | 129 | 122 | 116 | 110 | 104 | 97 | 91 |

---

## Usage Examples

```python
import numpy as np
from constants import alt_eff, nu_i, kappa_e, Psi_0

run = 5  # 104 km
print(f"h = {alt_eff[run]} km")
print(f"nu_i = {nu_i[run]:.3g} s^-1, kappa_e = {kappa_e[run]:.3g}, Psi_0 = {Psi_0[run]:.3g}")
```

Find the run index for a given altitude:

```python
import numpy as np
from constants import alt_eff

run = int(np.where(alt_eff == 110)[0][0])  # -> 4
```

!!! tip "Cache invalidation"
    [`grid_cache.py`](grid_cache.md) tracks the modification time of `constants.py`: if you edit any constant, previously cached CSV grids are automatically treated as stale and recomputed on the next run.
