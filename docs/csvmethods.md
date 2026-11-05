# `csvmethods.py` — Post-Processing Cached CSVs

## Overview

Analysis utilities that operate **directly on the cached long-format CSVs** in the repository's `data/` directory (as written by [`grid_cache.py`](grid_cache.md)), without re-solving anything. They support two physics questions:

1. **Where does the wave transition between regimes?** `plot_ratio_contour` plots $\mathrm{Im}(\omega)/\mathrm{Re}(\omega)$ but *only* inside the resonance band $|\omega| \approx \nu_{in}$ (the ion-neutral collision frequency), i.e. where the wave frequency is comparable to the collision rate — the demagnetized↔collisional transition region. The relevant $\nu_{in}$ is looked up automatically from the altitude encoded in the filename (e.g. `110.csv` → 110 km → `run` with `alt_eff == 110`).
2. **How much does the ion acceleration term matter?** `plot_field_difference` takes the pointwise difference of a field (default: growth rate) between two CSVs sampled on the same grid — e.g. an `eq_73` run vs. an `eq_73_ai0zero` run — producing a map of $\Delta\,\mathrm{Im}(\omega)$ attributable to $\mathbf{a}_{i0}$.

CSV filenames listed in `FILES` and processed in `__main__` follow the convention `data/<altitude>.csv` / `data/<altitude>_ai0zero.csv` (create them with `main.py --csv-name`).

---

## API Reference

### `load_ratio_grid(path)`

Pivot a long-format CSV into 2-D grids.

| Argument | Type | Description |
|---|---|---|
| `path` | `str` | CSV with columns `kx_m^-1`, `ky_m^-1`, `growth_rate_s^-1`, `frequency_rad_s^-1`. |

**Returns:** `(kx_vals, ky_vals, ratio_grid, omega_mag_grid)` — 1-D axes plus 2-D arrays of $\mathrm{Im}(\omega)/\mathrm{Re}(\omega)$ (±inf replaced by `nan`) and $|\omega| = \sqrt{\mathrm{Re}^2 + \mathrm{Im}^2}$.

### `max_ratio(path)`

| Argument | Type | Description |
|---|---|---|
| `path` | `str` | CSV path. |

**Returns:** `float` — maximum of $|\mathrm{Im}(\omega)/\mathrm{Re}(\omega)|$ over the grid (nan-aware).

### `nu_in_for_altitude(path)`

Look up $\nu_i$ for the altitude encoded as the first 3-digit number in the filename (e.g. `110_smallk.csv` → 110 km).

**Returns:** `float` — ion-neutral collision frequency (s⁻¹) for that altitude's run.

### `plot_ratio_contour(path, save_path, band_frac=0.2)`

Plot $\mathrm{Im}(\omega)/\mathrm{Re}(\omega)$ restricted to the resonance band $|\omega| \in [\,\nu_{in}(1-\text{band\_frac}),\ \nu_{in}(1+\text{band\_frac})\,]$, with a symmetric diverging color scale clipped at the 95th percentile of $|{\rm ratio}|$ and the band boundary contours hatched/labeled. Prints band statistics (point count, median and max $|\mathrm{Im}/\mathrm{Re}|$).

| Argument | Type | Description |
|---|---|---|
| `path` | `str` | Input CSV (altitude encoded in filename). |
| `save_path` | `str` | Output figure path. |
| `band_frac` | `float` | Half-width of the accepted $\lvert\omega\rvert \approx \nu_{in}$ band as a fraction of $\nu_{in}$. Default `0.2` (±20 %). |

**Returns:** `None` (saves a figure, or prints a message and returns if no points fall in the band).

### `load_field_grid(path, column="growth_rate_s^-1")`

Generic pivot of one CSV column into a 2-D grid.

**Returns:** `(kx_vals, ky_vals, grid)`.

### `plot_field_difference(path_nonzero, path_zero, save_path, column="growth_rate_s^-1", percentile=95)`

Contour of `path_nonzero` − `path_zero` for the chosen column, with a symmetric color scale clipped at `percentile` of $|{\rm diff}|$ so a few resonance outliers don't wash out the structure.

| Argument | Type | Description |
|---|---|---|
| `path_nonzero` | `str` | CSV of the full run (e.g. $a_{i0} \ne 0$). |
| `path_zero` | `str` | CSV of the reference run (e.g. $a_{i0} = 0$). Must share the exact same $k$ grid. |
| `save_path` | `str` | Output figure path. |
| `column` | `str` | CSV column to difference. Default growth rate. |
| `percentile` | `float` | Color-scale clip percentile of $\lvert\text{diff}\rvert$. Default `95`. |

**Returns:** `None`. **Raises** `ValueError` if the grids differ or no finite differences exist.

---

## Usage Examples

Produce the input CSVs with friendly names, then post-process:

```bash
# 1. Generate matched a_i0-on / a_i0-off grids at 104 km
python main.py --method growth_contour --save tmp_on.png  --run 5 --csv-name 104.csv
python main.py --equation eq_73_ai0zero --method growth_contour \
    --save tmp_off.png --run 5 --csv-name 104_ai0zero.csv

# 2. Run the built-in difference plots (data/104*.csv, data/110*.csv)
python csvmethods.py
```

From Python:

```python
from csvmethods import plot_ratio_contour, plot_field_difference, max_ratio

print("max |Im/Re| at 104 km:", max_ratio("data/104.csv"))

# Im/Re ratio inside the |omega| ~ nu_in transition band (±30%)
plot_ratio_contour("data/104.csv", "plots/ratio_104.png", band_frac=0.3)

# Effect of the ion acceleration term on the growth rate
plot_field_difference("data/104.csv", "data/104_ai0zero.csv", "plots/difference_104.png")
```
