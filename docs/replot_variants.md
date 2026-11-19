# `replot_variants.py` — Batch Regeneration of Variant Grids & Plots

## Overview

Script that regenerates the hand-named CSV grids in `data/` and the corresponding figures in `plots/` for every combination of altitude and physics variant:

| Variant | Equation | E-field scale | CSV | Plot |
|---|---|---|---|---|
| `normal` | `dispersion_rel` | ×1 | `data/<alt>.csv` | `plots/<alt>.png` |
| `ai0zero` | `dispersion_rel_ai0zero` | ×1 | `data/<alt>_ai0zero.csv` | `plots/<alt>_ai0zero.png` |
| `doubleE` | `dispersion_rel` | ×2 | `data/<alt>_doubleE.csv` | `plots/<alt>_double_E.png` |
| `doubleE_ai0zero` | `dispersion_rel_ai0zero` | ×2 | `data/<alt>_doubleE_ai0zero.csv` | `plots/<alt>_double_E_ai0zero.png` |
| `reduced` | `dispersion_rel_reduced` | ×1 | `data/<alt>_reduced.csv` | `plots/<alt>_reduced.png` |
| `reduced_ai0zero` | `dispersion_rel_reduced` | ×1 | `data/<alt>_reduced_ai0zero.csv` | `plots/<alt>_reduced_ai0zero.png` |
| `reduced_doubleE` | `dispersion_rel_reduced` | ×2 | `data/<alt>_reduced_doubleE.csv` | `plots/<alt>_reduced_double_E.png` |

The `reduced*` variants use [`dispersion_rel_reduced.py`](dispersion_rel_reduced.md) (fluid electrons, ion term with the $i\,\mathbf{k}\cdot\mathbf{a}_{i0}$ terms dropped). Since that equation no longer contains $\mathbf{a}_{i0}$, `reduced_ai0zero` solves the same equation as `reduced` and exists only for filename parity with the full set.

`<alt>` is the altitude in km: run 4 → `110`, run 5 → `104` (the `RUNS` map).

Design points:

- **One subprocess per (run, variant).** Doubling the electric field works by setting the `E_SCALE` environment variable *before* [`constants.py`](constants.md) is imported, so that every derived constant ($T_e$, $\lambda_D$, drifts, $\mathbf{a}_{i0}$, …) is consistently rescaled — exactly as if `V_d` had been edited by hand. Since `constants.py` computes these at import time, each combination must run in a fresh Python process.
- **Positive-only plots.** Whenever the chosen plotting method supports it (e.g. the default `growth_contour`), the script passes `positive_only=True` to [`main.run`](main.md), so only unstable regions ($\gamma > 0$) are colored. For methods without `positive_only` support (e.g. `growth_damping`) it prints a note and plots all values.
- **Always recomputes** (`force=True`) — this script exists to refresh the caches after physics changes.

---

## API Reference

### `replot_all_variants(method="growth_contour", runs=None, variants=None, points=100, k_max=100.0)`

Spawn one worker subprocess per (run, variant) combination. `runs` / `variants` default to all keys of `RUNS` / `VARIANTS`.

### CLI

| Flag | Default | Description |
|---|---|---|
| `--method NAME` | `growth_contour` | Plotting method from [`methods.py`](methods.md). |
| `--runs N [N ...]` | `4 5` | alt_eff indices (4 = 110 km, 5 = 104 km). |
| `--variants NAME [NAME ...]` | all seven | Subset of `VARIANTS`. |
| `--points N` | `100` | Grid resolution per axis. |
| `--k-max X` | `100.0` | Grid spans ±`k_max` in $k_x$ and $k_y$ (m⁻¹). |

(`--worker` is internal: it marks the per-combination subprocess.)

---

## Usage Examples

```bash
python replot_variants.py                          # everything, both altitudes
python replot_variants.py --variants doubleE       # just the doubled field
python replot_variants.py --points 20              # coarse grid for a quick look
python replot_variants.py --method growth_damping  # two-panel Im/Re plots instead
python replot_variants.py --variants reduced reduced_ai0zero reduced_doubleE  # reduced equation only
```
