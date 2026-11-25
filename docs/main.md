# `main.py` — Command-Line Entry Point

## Overview

`main.py` is the orchestrator: it wires a chosen **dispersion equation** ([`dispersion_rel.py`](dispersion_rel.md), [`dispersion_rel_reduced.py`](dispersion_rel_reduced.md) or [`electron_kinetic.py`](electron_kinetic.md)) to a chosen **plotting method** ([`methods.py`](methods.md)). Running it:

1. builds linear $k_x$ / $k_y$ grids from the CLI bounds,
2. solves the growth/frequency grid via [`root.solve_grid`](root.md),
3. hands the grid to the selected plotting method, which saves a figure.

Two registries map CLI names to callables:

| Registry | Keys |
|---|---|
| `EQUATIONS` | `dispersion_rel` (default), `dispersion_rel_ai0zero` — fluid electrons ([`dispersion_rel.py`](dispersion_rel.md)); `dispersion_rel_reduced` — fluid electrons, ion term without the $i\,\mathbf{k}\cdot\mathbf{a}_{i0}$ terms ([`dispersion_rel_reduced.py`](dispersion_rel_reduced.md)); `eq_72_kinetic`, `eq_72_kinetic_ae0zero` (electron NTA limit, $a_{e0}=0$), `eq_72_kinetic_ai0zero` — fully kinetic electrons ([`electron_kinetic.py`](electron_kinetic.md)) |
| `METHODS` (from `methods.py`) | `growth_damping`, `growth_contour`, `growth_frequency_ratio`, `phase_velocity` |

---

## API Reference

### `run(equation, method, save_path, run=5, kx_min=-100.0, kx_max=100.0, ky_min=-100.0, ky_max=100.0, points=100, positive_only=False)`

Programmatic equivalent of the CLI.

| Argument | Type | Description |
|---|---|---|
| `equation` | `str` or `callable` | Key in `EQUATIONS`, or a dispersion function directly. |
| `method` | `str` | Key in `METHODS`. |
| `save_path` | `str` | Output figure filename. |
| `run` | `int` | Altitude index into `alt_eff`. Default `5` (104 km). |
| `kx_min`, `kx_max`, `ky_min`, `ky_max` | `float` | Wave-number bounds (m⁻¹). Default ±100. |
| `points` | `int` | Grid points per axis. Default `100`. |
| `positive_only` | `bool` | Mask non-positive growth rates (only supported by `growth_contour`; raises `ValueError` otherwise). |

**Returns:** `None` — the side effect is a saved figure.

### `main()`

Parses the command-line arguments below and calls `run`.

| Flag | Default | Description |
|---|---|---|
| `--equation NAME` | `dispersion_rel` | Dispersion relation to solve — any key of `EQUATIONS` (see Overview). |
| `--method {growth_contour, growth_damping, growth_frequency_ratio, phase_velocity}` | *(required)* | Plotting method. |
| `--save PATH` | *(required)* | Output plot filename. |
| `--run N` | `5` | Index into `alt_eff` (0 → 135 km … 7 → 91 km). |
| `--kx-min/--kx-max/--ky-min/--ky-max X` | ±100.0 | Wave-number bounds (m⁻¹). |
| `--points N` | `100` | Grid resolution per axis. |
| `--positive-only` | off | Mask non-positive growth rates (`growth_contour` only). |

---

## Usage Examples

```bash
# Two-panel Im/Re(omega) plot at 104 km (run 5)
python main.py --method growth_damping --save growth_damping_104km.png --run 5

# Growth-rate contour at 116 km, only unstable regions
python main.py --method growth_contour --save growth_116km.png --run 3 --positive-only

# Higher resolution, smaller k window
python main.py --method phase_velocity --save vph_110km.png --run 4 \
    --kx-min -40 --kx-max 40 --ky-min -40 --ky-max 40 --points 200

# a_i0 = 0 comparison run
python main.py --equation dispersion_rel_ai0zero --method growth_contour \
    --save growth_104km_ai0zero.png --run 5

# Fully kinetic electrons (manuscript Eq. 72/75)
python main.py --equation eq_72_kinetic --method growth_contour \
    --save growth_104km_kinetic.png --run 5

# Reduced equation (ion term without the i k.a_i0 terms)
python main.py --equation dispersion_rel_reduced --method growth_contour \
    --save growth_104km_reduced.png --run 5
```

From Python:

```python
import main

main.run("dispersion_rel", "growth_contour", "growth_contour_116km.png", run=3)
```
