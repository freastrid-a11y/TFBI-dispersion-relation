# `root.py` — Root Finding & Grid Solving

## Overview

This module turns the dispersion residual from [`eq73.py`](eq73.md) into actual solutions $\omega(\mathbf{k})$. It provides:

1. **`Z_tilde`** — the modified plasma dispersion function of manuscript Eq. (74), built on the Faddeeva function $w(z)$ (`scipy.special.wofz`). The branch is selected by the sign of $\mathrm{Im}(\omega_i)$, where $\omega_i = \omega + i\nu_i$ is the collision-shifted frequency, so that the analytic continuation of the velocity-space integral is taken consistently for **both** ion arguments $z_1$ and $z_2$.

2. **A robust multistart root finder** (`find_roots`). The dispersion relation is stiff, has many spurious near-flat regions, and its natural frequency scale varies by orders of magnitude across the $(k_x, k_y)$ plane. The solver therefore:
   - nondimensionalizes $\omega$ by a characteristic `frequency_scale` (max of $\nu_i$, $k|V_d|$, $kV_{Tx}$, $kV_{Ty}$);
   - normalizes the residual by a fixed electron-reference scale so one tolerance is meaningful at every $k$;
   - launches `scipy.optimize.root` (Powell hybrid, `method="hybr"`) from a grid of physically motivated initial guesses (drift-wave frequency, electron/ion drift frequencies, ± fractions of the scale);
   - **independently verifies** each candidate: finite residual below `residual_tol`, estimated remaining Newton correction below `correction_tol`, and $|\omega|$ within `max_frequency_factor` of the scale (rejects spurious convergence at huge $|\omega|$);
   - merges duplicates reached from different seeds.

3. **A continuation-based grid sweep** (`solve_grid`) that marches through the $(k_x, k_y)$ grid, seeding each point with the converged roots of its already-solved neighbors so the selected branch varies smoothly.

With the convention $\exp(i\mathbf{k}\cdot\mathbf{r} - i\omega t)$, $\mathrm{Im}(\omega) > 0$ is **growth** and $\mathrm{Re}(\omega)$ is the real wave frequency.

---

## API Reference

### `Z_tilde(z, omega_i)`

Modified plasma dispersion function (manuscript Eq. 74).

| Argument | Type | Description |
|---|---|---|
| `z` | `complex` | Argument of the dispersion function ($z_1$ or $z_2$ in `eq_73`). |
| `omega_i` | `complex` | Collision-shifted frequency $\omega + i\nu_i$; its imaginary sign selects the branch. |

**Returns:** `complex` —
$i\sqrt{\pi}\,w(z)$ for $\mathrm{Im}(\omega_i) > 0$; $-i\sqrt{\pi}\,w(-z)$ for $\mathrm{Im}(\omega_i) < 0$; `nan+nan*j` exactly on the real axis (the manuscript leaves this case unspecified, and the solver discards non-finite evaluations).

### `frequency_scale(k_x, k_y, run)`

Characteristic frequency used to nondimensionalize $\omega$ for the solver.

| Argument | Type | Description |
|---|---|---|
| `k_x`, `k_y` | `float` | Wave-vector components (m⁻¹). |
| `run` | `int` | Altitude index into `constants.py` arrays. |

**Returns:** `float` — $\max\!\left(1,\ \nu_i,\ k|V_d|,\ kV_{Tx},\ kV_{Ty}\right)$ in s⁻¹.

### `default_guesses(k_x, k_y, run, scale)`

Build the multistart seed list: the cross product of several real-frequency guesses (0, the combined drift-wave estimate $\omega_\mathrm{drift} = (\omega_e + \Psi_0\,\omega_i)/(1+\Psi_0)$, pure electron/ion drift frequencies, $\pm 0.5\,$scale) with three growth-rate guesses ($+0.05$, $-0.05$, $-0.5$ × scale).

| Argument | Type | Description |
|---|---|---|
| `k_x`, `k_y` | `float` | Wave-vector components (m⁻¹). |
| `run` | `int` | Altitude index. |
| `scale` | `float` | Output of `frequency_scale`. |

**Returns:** `list[complex]` — initial guesses for $\omega$.

### `find_roots(dispersion, k_x, k_y, run, guesses=None, include_default_guesses=True, residual_tol=1e-8, correction_tol=1e-6, max_frequency_factor=30.0)`

Find accepted roots of `dispersion(omega, k_x, k_y, run) = 0` at one wave vector.

| Argument | Type | Description |
|---|---|---|
| `dispersion` | `callable` | Dispersion residual, e.g. `eq_73`. Signature `(omega, k_x, k_y, run) -> complex`. |
| `k_x`, `k_y` | `float` | Wave-vector components (m⁻¹). Returns `[]` at exactly $k = 0$. |
| `run` | `int` | Altitude index. |
| `guesses` | `list[complex]` or `None` | Extra seeds (e.g. neighboring grid solutions for continuation). |
| `include_default_guesses` | `bool` | Also use `default_guesses`. Default `True`. |
| `residual_tol` | `float` | Max accepted normalized residual \|ε\|. Default `1e-8`. |
| `correction_tol` | `float` | Max accepted estimated Newton correction, relative to the frequency scale. Default `1e-6`. |
| `max_frequency_factor` | `float` | Reject roots with $\|\omega\| >$ factor × scale. Default `30.0`. |

**Returns:** `list[complex]` — distinct accepted roots, **sorted by decreasing $\mathrm{Im}(\omega)$** (most unstable first).

!!! note
    This is a *local, multistart* solver — it is not a proof that all roots were found.

### `solve_grid(dispersion, kx_vals, ky_vals, run=0)`

Sweep a wave-vector grid, keeping the root with the largest $\mathrm{Im}(\omega)$ at each point and seeding each point from its left and lower neighbors (branch continuation). Prints one progress line per row.

| Argument | Type | Description |
|---|---|---|
| `dispersion` | `callable` | Dispersion residual, e.g. `eq_73`. |
| `kx_vals`, `ky_vals` | `array_like` of `float` | 1-D grids of wave numbers (m⁻¹). |
| `run` | `int` | Altitude index. Default `0`. |

**Returns:** `ndarray[complex]` of shape `(len(ky_vals), len(kx_vals))` — the dominant root $\omega$ at each grid point; `nan+nan*j` where no root was accepted.

---

## Usage Examples

Find the dominant mode at one wave vector:

```python
from eq73 import eq_73
from root import find_roots

roots = find_roots(eq_73, k_x=40.0, k_y=0.0, run=5)
if roots:
    omega = roots[0]
    print(f"Re(omega) = {omega.real:.3g} rad/s, growth = {omega.imag:.3g} 1/s")
```

Solve a small grid directly (bypassing the CSV cache):

```python
import numpy as np
from eq73 import eq_73
from root import solve_grid

kx = np.linspace(-50, 50, 25)
ky = np.linspace(-50, 50, 25)
omega_grid = solve_grid(eq_73, kx, ky, run=5)

growth = omega_grid.imag      # Im(omega) > 0 => unstable
frequency = omega_grid.real
```

In normal use you should go through the cache instead:

```bash
python main.py --method growth_damping --save growth_damping_104km.png --run 5
```
