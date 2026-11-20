# `electron_kinetic.py` — Fully Kinetic Electrons (Eqs. 67, 70, 72, 75)

## Overview

This module replaces the fluid electron susceptibility of [`dispersion_rel.py`](dispersion_rel.md) with the **fully kinetic electron contribution** of the manuscript — Eqs. (67a–c), (70) and (75) — and assembles the full kinetic dispersion relation, Eq. (72):

$$
I_\mathrm{ion} + I_\mathrm{electron} + k^2\lambda_D^2 \;=\; 0 .
$$

The result has the **same overall normalisation** as `dispersion_rel` (which returns $k^2\lambda_D^2(1+\chi_e+\chi_i)$), so the root finder in [`root.py`](root.md) and the pipeline in [`main.py`](main.md) work unchanged.

!!! warning "Notation: the paper's electron $k_x$, $k_y$ are not this repo's $k_x$, $k_y$"
    The manuscript uses two different meanings of $k_x, k_y$:

    - **Ion part, Eq. (73):** $x$ is along $\mathbf{a}_{i0}$ (i.e. $E_0$), and $\mathbf{k}$ lies in the plane perpendicular to $\mathbf{B}$.
    - **Electron part, Eqs. (23)–(25), (43)–(46), (67), (75):** the paper's "$k_y$" is $k_\perp$ (the magnitude of $\mathbf{k}$ perpendicular to $\mathbf{B}$), while "$k_x$" is $k_\parallel$ (parallel to $\mathbf{B}$).

    This repository always uses the frame: $\mathbf{V}_d$ ($\mathbf{E}\times\mathbf{B}$ drift) along $+x$, $E_0$ (and $\mathbf{a}_{i0}$, $\mathbf{a}_{e0}$) along $+y$, $\mathbf{B}$ along $+z$. In this module, therefore:

    | Paper | This module |
    |---|---|
    | $k_y^2$ (electron part) | $k_\perp^2 = k_x^2 + k_y^2$ |
    | $k_x$ (electron part) | `k_par` (0 for strictly perpendicular waves) |
    | $\mathbf{k}\cdot\mathbf{a}_{e0}$ | $k_y\, a_{e0}$ |

Other definitions used here (paper Eqs. 29–30): $a_{e0} = eE_0/m_e$, $V_{Te}^2 = T_e/m_e$, $\nu_0 = \langle\nu_{en}\rangle/3$ (since $\nu_{en}(V) = \nu_0 V^2/V_{Te}^2$), $\Delta_{\omega k} = \mathbf{k}\cdot\mathbf{V}_{e0} - \omega$, and the ordinary plasma dispersion function $Z(x) = i\sqrt{\pi}\,w(x)$.

---

## API Reference

### `Z(x)`

Ordinary plasma dispersion function $Z(x) = i\sqrt{\pi}\,w(x)$ (via `scipy.special.wofz`).

### `first_quadrant_quartic_root(q4)`

Return the 4th root of `q4` lying in the open first quadrant ($\mathrm{Re}\,q > 0$, $\mathrm{Im}\,q > 0$) — the $q_0$ of Eq. (67c). The four candidates are $q_4^{1/4} i^n$; the root that maximises $\min(\mathrm{Re}\,q, \mathrm{Im}\,q)$ is chosen, which is unique and continuous. Accepts scalars or arrays.

### `electron_params(omega, k_x, k_y, run, k_par=0.0, a_e0_enabled=True)`

Compute the building blocks of Eqs. (67a–c).

| Argument | Type | Description |
|---|---|---|
| `omega` | `complex` | Complex wave frequency (rad/s). |
| `k_x`, `k_y` | `float` | Wave numbers in this repo's frame (m⁻¹). |
| `run` | `int` | Altitude index into [`constants.py`](constants.md) arrays. |
| `k_par` | `float` | Field-aligned wave number (paper's electron "$k_x$"). Default `0`. |
| `a_e0_enabled` | `bool` | If `False`, sets $a_{e0} = 0$ (NTA limit). Default `True`. |

**Returns:** `dict` with keys `D1`, `alpha0`, `alpha1`, `q4`, `nu_0`, `Delta`, `V_Te_sq`, `a_e0`, `k_dot_a`, where

$$
D_1 = 2\left(k_\perp^2 V_{Te}^2 - i\,\mathbf{k}\cdot\mathbf{a}_{e0}\right),
$$

$\alpha_1 = 5i\,\mathbf{k}\cdot\mathbf{a}_{e0}/(2 D_1)$ (67a), $\alpha_0 = k_\parallel^2 V_{Te}^2\,\Omega_e^2/(4\nu_0^2 D_1)$ (67b), and $q^4 = -(3i\nu_0\Delta_{\omega k} + k_\parallel^2 V_{Te}^2)\,\Omega_e^2/(4\nu_0^2 D_1)$ (67c).

### `I_electron(omega, k_x, k_y, run, k_par=0.0, a_e0_enabled=True)`

Kinetic electron contribution $I_\mathrm{electron}$, manuscript Eq. (75) — the bracket alone, with no $D_1/D_2$ prefactor:

$$
I_e = 1 + 2\alpha_1
+ \frac{q^4 + \alpha_0 + \alpha_1 q_0^2}{q_0}\,Z(q_0)
- i\,\frac{q^4 + \alpha_0 - \alpha_1 q_0^2}{q_0}\,Z(iq_0).
$$

**Returns:** `complex` (or array, broadcasting over the inputs).

### `I_ion(omega, k_x, k_y, run, a_i0_enabled=True)`

Ion contribution $I_\mathrm{ion}$, manuscript Eq. (73), in this repo's frame ($\mathbf{a}_{i0}$ and $\mathbf{U}_0$ along $+y$). Identical physics to `dispersion_rel`'s $\chi_i$, multiplied by $k^2\lambda_D^2$.

!!! note "Eq. (73) prefactor typo"
    The prefactor uses $\mathrm{den}_2$ (the bi-Maxwellian denominator), as in manuscript Eq. (15). Eq. (73) as printed has $\mathrm{den}_1$ there, which is a typo. The same correction is applied in [`dispersion_rel.py`](dispersion_rel.md).

### `eq_72_kinetic(omega, k_x, k_y, run, a_i0_enabled=True, a_e0_enabled=True, k_par=0.0)`

Full kinetic dispersion function, Eq. (72): $I_\mathrm{ion} + I_\mathrm{electron} + k^2\lambda_D^2$ ($k^2$ floored at $10^{-12}$). Registered in `main.EQUATIONS` as `eq_72_kinetic`.

### `eq_72_kinetic_ae0zero(omega, k_x, k_y, run)`

Eq. (72) in the **NTA limit for electrons** ($a_{e0} \to 0$), with $\mathbf{a}_{i0}$ kept. Registered as `eq_72_kinetic_ae0zero`.

### `eq_72_kinetic_ai0zero(omega, k_x, k_y, run)`

Eq. (72) with $\mathbf{a}_{i0} \to 0$, kinetic electrons with $a_{e0}$ kept. Registered as `eq_72_kinetic_ai0zero`.

---

## Usage Examples

```bash
# Growth-rate contour with fully kinetic electrons at 104 km
python main.py --equation eq_72_kinetic --method growth_contour \
    --save growth_104km_kinetic.png --run 5 --csv-name 104_kinetic.csv
```

From Python:

```python
from electron_kinetic import eq_72_kinetic
from root import find_roots

roots = find_roots(eq_72_kinetic, k_x=30.0, k_y=10.0, run=5)
print(roots[0])  # fastest-growing kinetic-electron mode
```

For a side-by-side fluid vs. kinetic comparison, see [`compare_electrons.py`](compare_electrons.md).
