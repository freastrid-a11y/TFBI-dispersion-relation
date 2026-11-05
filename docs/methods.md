# `methods.py` — Plotting Methods

## Overview

Visualization layer. Every public function here consumes an **already-solved** grid — as produced and cached by [`grid_cache.get_grid`](grid_cache.md) — and saves a figure. All methods share the signature

```python
method(kx_vals, ky_vals, growth, frequency, altitude, save_path, **options)
```

where `growth` = $\mathrm{Im}(\omega)$ and `frequency` = $\mathrm{Re}(\omega)$, both 2-D arrays of shape `(len(ky_vals), len(kx_vals))` with `nan` at points where no root was accepted (rendered light gray / masked).

The module exports the `METHODS` registry used by `main.py`:

| Key | Function | What it shows |
|---|---|---|
| `growth_damping` | `growth_damping` | Two panels: $\mathrm{Im}(\omega)$ and $\mathrm{Re}(\omega)$, diverging colormap centered on 0, with the marginal-stability contour $\gamma = 0$ overlaid. |
| `growth_contour` | `growth_contour` | Single-panel growth-rate map (viridis); optionally only $\gamma > 0$. |
| `growth_frequency_ratio` | `growth_frequency_ratio` | $\mathrm{Im}(\omega)/\mathrm{Re}(\omega)$ — a normalized measure of how strongly a wave grows per oscillation. |
| `phase_velocity` | `phase_velocity` | $\mathrm{Re}(\omega)/\lvert k\rvert$ — the wave phase speed (m/s), comparable to the ion-acoustic speed $C_s$ in FBI theory. |

Physics note: the black $\gamma = 0$ contour in `growth_damping` traces the **marginal stability boundary** — inside it the TFBI grows, outside it the modes are damped.

---

## API Reference

### `growth_damping(kx_vals, ky_vals, growth, frequency, altitude, save_path)`

Two-panel $\mathrm{Im}(\omega)$ / $\mathrm{Re}(\omega)$ pcolormesh with a shared `RdBu_r` diverging colormap centered on zero (`TwoSlopeNorm`), symmetric tick layout, invalid points in light gray, and the $\gamma = 0$ contour on both panels.

| Argument | Type | Description |
|---|---|---|
| `kx_vals`, `ky_vals` | `ndarray[float]` | 1-D wave-number grids (m⁻¹). |
| `growth` | `ndarray[float]` 2-D | $\mathrm{Im}(\omega)$ (s⁻¹). |
| `frequency` | `ndarray[float]` 2-D | $\mathrm{Re}(\omega)$ (s⁻¹). |
| `altitude` | `float` | Effective altitude (km) for the title. |
| `save_path` | `str` | Output figure path (saved at 250 dpi). |

**Returns:** `None` (saves and closes the figure, prints the path).

### `growth_contour(kx_vals, ky_vals, growth, frequency, altitude, save_path, positive_only=False)`

Single-panel growth-rate contour (viridis).

Additional argument:

| Argument | Type | Description |
|---|---|---|
| `positive_only` | `bool` | If `True`, masks all $\gamma \le 0$ so only unstable regions are colored. Exposed via `main.py --positive-only`. |

**Returns:** `None`. If the mask leaves no data, an explanatory "No accepted (unstable) roots found" message is drawn instead.

### `growth_frequency_ratio(kx_vals, ky_vals, growth, frequency, altitude, save_path)`

$\mathrm{Im}(\omega)/\mathrm{Re}(\omega)$ contour (`RdBu_r`), undefined (masked) where $\mathrm{Re}(\omega) = 0$ or either value is non-finite.

**Returns:** `None`.

### `phase_velocity(kx_vals, ky_vals, growth, frequency, altitude, save_path)`

Phase-velocity contour $\mathrm{Re}(\omega)/\sqrt{k_x^2 + k_y^2}$ in m/s (`RdBu_r`), undefined at $k = 0$.

**Returns:** `None`.

### Internal helpers

- `_actual_range_norm(data)` — builds a `TwoSlopeNorm` spanning the actual finite data range with center 0 (guards against degenerate all-positive/all-negative data).
- `_add_zero_growth_boundary(ax, kx_vals, ky_vals, growth)` — overlays and labels the $\gamma = 0$ contour wherever the data straddles zero.

---

## Usage Examples

Via the CLI (recommended — uses the CSV cache):

```bash
python main.py --method growth_damping --save growth_damping_104km.png --run 5
python main.py --method growth_contour --save growth_104km.png --run 5 --positive-only
python main.py --method phase_velocity --save vph_104km.png --run 5
```

Direct use on a grid you already have:

```python
import numpy as np
from eq73 import eq_73
from grid_cache import get_grid
from methods import growth_frequency_ratio
from constants import alt_eff

run = 4  # 110 km
kx = np.linspace(-100, 100, 100)
ky = np.linspace(-100, 100, 100)
kx, ky, growth, frequency = get_grid(eq_73, kx, ky, run)

growth_frequency_ratio(kx, ky, growth, frequency, alt_eff[run], "ratio_110km.png")
```
