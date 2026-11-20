# TFBI Dispersion Relation Solver

Numerical solver for the **Thermal Farley–Buneman Instability (TFBI)** dispersion relation in the lower ionosphere (E region). It finds complex roots $\omega(\mathbf{k})$ over a 2-D $(k_x, k_y)$ grid at a chosen effective altitude, caches the results as CSV, and plots growth rate $\mathrm{Im}(\omega)$, frequency $\mathrm{Re}(\omega)$ and related quantities.

The equations come from the kinetic TFBI theory manuscript (Dimant et al.), which is in `papers/`. The solver includes:

- **Fluid electrons + kinetic ions** (manuscript Eq. 73), with and without ion acceleration $\mathbf{a}_{i0}$
- **Reduced** form: fluid electrons, with the $i\,\mathbf{k}\cdot\mathbf{a}_{i0}$ terms dropped from the ion term
- **Fully kinetic electrons + kinetic ions** (Eqs. 67, 70, 72, 75), with optional $\mathbf{a}_{i0} = 0$ / $\mathbf{a}_{e0} = 0$ limits

📖 **Full documentation:** https://freastrid-a11y.github.io/TFBI-dispersion-relation/ (built from `docs/` with MkDocs)

---

## Requirements

Python 3.10+ with:

```bash
pip install numpy scipy pandas matplotlib
```

To build the docs locally you also need `pip install mkdocs-material`.

## Quick start

Solve the fluid-electron dispersion relation at 104 km (`--run 5`) on a 100×100 grid and plot growth and frequency:

```bash
python main.py --method growth_damping --save growth_damping_104km.png --run 5
```

Use fully kinetic electrons instead:

```bash
python main.py --equation eq_72_kinetic --method growth_contour --save 104_kinetic.png --run 5
```

The first run solves the grid and writes it to `data/`. Later runs with the same setup read the cached CSV. The cache is recomputed automatically if `constants.py` or the equation's source file has changed since it was written; pass `--force` to recompute anyway.

### `main.py` options

| Flag | Default | Description |
|---|---|---|
| `--equation` | `dispersion_rel` | One of `dispersion_rel`, `dispersion_rel_ai0zero`, `dispersion_rel_reduced`, `eq_72_kinetic`, `eq_72_kinetic_ai0zero`, `eq_72_kinetic_ae0zero` |
| `--method` | *(required)* | `growth_damping`, `growth_contour`, `growth_frequency_ratio`, `phase_velocity` |
| `--save` | *(required)* | Output image path |
| `--run` | `5` | Index into `alt_eff = [135, 129, 122, 116, 110, 104, 97, 91]` km (4 → 110 km, 5 → 104 km) |
| `--kx-min/--kx-max`, `--ky-min/--ky-max` | `±100` | Wavevector range (m⁻¹) |
| `--points` | `100` | Grid points per axis |
| `--positive-only` | off | Mask non-positive growth (`growth_contour` only) |
| `--force` | off | Ignore the cache and re-solve |
| `--csv-name` | auto | Name for the cached CSV in `data/` |

## Batch scripts

```bash
# Regenerate every variant grid + plot (normal, a_i0 = 0, doubled E, reduced, kinetic) for 104 and 110 km
python replot_variants.py
python replot_variants.py --variants doubleE --runs 4      # one variant, one altitude
python replot_variants.py --points 20                      # coarse grid for a quick look

# Four-panel comparison: reduced (old) vs kinetic Eq. 72 (new), at E_0 and 2E_0
python fourpanel.py --alt 110
python fourpanel.py --alt 104 --separate-scales
```

The doubled-field variants set the `E_SCALE` environment variable, which `constants.py` reads, so every derived quantity ($T_e$, $\lambda_D$, drifts, $a_{i0}$, …) is rescaled consistently.

`csvmethods.py` post-processes cached CSVs, for example resonance bands in $\mathrm{Im}\,\omega/\mathrm{Re}\,\omega$ and growth-rate differences between $a_{i0} \neq 0$ and $a_{i0} = 0$.

## Repository layout

| Path | Contents |
|---|---|
| `constants.py` | Physical constants, background fields, altitude-dependent collision profile ($\psi(h)$, $\kappa_{i,e}$, $\nu_{i,e}$) |
| `dispersion_rel.py` | Fluid-electron / kinetic-ion dispersion relation |
| `dispersion_rel_reduced.py` | Reduced form with the $i\,\mathbf{k}\cdot\mathbf{a}_{i0}$ terms dropped |
| `electron_kinetic.py` | Kinetic electron and ion terms and the full kinetic relation (Eq. 72) |
| `root.py` | Modified plasma dispersion function $\tilde Z$, multistart Newton root finder with neighbor passing, grid sweep |
| `grid_cache.py` | Solve-once-and-cache layer for grids in `data/` |
| `methods.py` | Plotting methods used by `main.py` |
| `main.py` | Command-line entry point |
| `replot_variants.py`, `fourpanel.py`, `csvmethods.py` | Batch regeneration, comparison figures, CSV post-processing |
| `claude_solution/` | Independent implementation of Eq. 72 (argument-principle root counting + Newton), used to cross-check the main solver (`solve.py`, `plot_indep.py`) |
| `data/` | Cached solution grids (CSV) |
| `plots/`, `compare/` | Generated figures and comparison outputs |
| `docs/`, `mkdocs.yml` | Per-module documentation, deployed to GitHub Pages on push to `main` |
| `papers/` | Reference manuscripts |
| `old_stuff/` | Earlier equation versions and their outputs, kept for reference |

## Cross-checking with the independent solver

```bash
cd claude_solution
python solve.py --alt 104                  # writes out/indep_104.npz
python plot_indep.py --alt 104 --compare   # compares with data/104_kinetic.csv
```

## Documentation

Serve the docs locally:

```bash
mkdocs serve
```
