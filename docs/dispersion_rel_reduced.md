# `dispersion_rel_reduced.py` — Reduced TFBI Dispersion Relation

## Overview

A reduced variant of the core relation in [`dispersion_rel.py`](dispersion_rel.md): the **fluid electron susceptibility $\chi_e$ is kept verbatim**, while the kinetic ion term is reduced by **dropping the $i\,\mathbf{k}\cdot\mathbf{a}_{i0}$ terms** from both kinetic denominators:

$$
\mathrm{den}_1 = k^2 V_{Tn}^2,
\qquad
\mathrm{den}_2 = k_x^2 V_{Tx}^2 + k_y^2 V_{Ty}^2 .
$$

Everything else — frictional ion heating, temperature anisotropy, the Pedersen drift Doppler shift $\mathbf{k}\cdot\mathbf{U}_0$, the $\tilde Z(z_1)$ collision correction, the den$_2$ prefactor (Eq. 15), and the $k^2\lambda_D^2$ numerical scaling — is identical to `dispersion_rel`.

!!! note "Equivalence to the $a_{i0} = 0$ limit"
    In the fluid-electron / kinetic-ion relation, $\mathbf{a}_{i0}$ enters **only** through the $i\,\mathbf{k}\cdot\mathbf{a}_{i0}$ denominators, so dropping them is mathematically identical to the NTA limit $\mathbf{a}_{i0} \to 0$: `dispersion_rel_reduced` returns exactly the same values as [`dispersion_rel_ai0zero`](dispersion_rel.md) (verified numerically). It is kept as a separate named equation so the reduced model gets its **own CSV caches and plot files** (`*_reduced*`), independent of the `ai0zero` comparison runs.

---

## API Reference

### `dispersion_rel_reduced(omega, k_x, k_y, run)`

Evaluate the scaled reduced dispersion residual $k^2\lambda_D^2\,(1 + \chi_e + \chi_i^{\mathrm{red}})$.

| Argument | Type | Description |
|---|---|---|
| `omega` | `complex` | Complex wave frequency $\omega = \omega_r + i\gamma$ (rad/s). |
| `k_x` | `float` | Wave number along the electron drift axis (m⁻¹). |
| `k_y` | `float` | Wave number along the electric-field / Pedersen axis (m⁻¹). |
| `run` | `int` | Index into the altitude-dependent arrays in [`constants.py`](constants.md). |

**Returns:** `complex` — the scaled dispersion residual; (numerically) zero when $(\omega, \mathbf{k})$ satisfies the reduced relation.

The interface matches `dispersion_rel`, so [`root.find_roots`](root.md) and [`main.run`](main.md) work unchanged. It is registered in `main.EQUATIONS` as `dispersion_rel_reduced`.

---

## Usage Examples

```bash
# Growth-rate contour of the reduced relation at 110 km
python main.py --equation dispersion_rel_reduced --method growth_contour \
    --save growth_110km_reduced.png --run 4
```

Batch-generate all reduced variant plots (see [`replot_variants.py`](replot_variants.md)):

```bash
python replot_variants.py --variants reduced reduced_ai0zero reduced_doubleE
```

From Python:

```python
from dispersion_rel_reduced import dispersion_rel_reduced
from root import find_roots

roots = find_roots(dispersion_rel_reduced, k_x=30.0, k_y=10.0, run=4)
print(roots)  # sorted by decreasing Im(omega)
```
