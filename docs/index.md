# TFBI Solver Documentation

Welcome to the documentation for the Thermal Farley-Buneman Instability (TFBI) solver.

## Getting Started
To run the solver, use the main entry point:

```bash
python main.py --method growth_damping --save output.png --run 5
```

## Module Guide

| Module | Role |
|---|---|
| [dispersion_rel](dispersion_rel.md) | The TFBI dispersion relation with fluid electrons (manuscript Eq. 73) — the core math. |
| [electron_kinetic](electron_kinetic.md) | Fully kinetic electron contribution and the full kinetic dispersion relation (Eqs. 67, 70, 72, 75). |
| [root](root.md) | Modified plasma dispersion function, multistart root finder, grid sweep. |
| [constants](constants.md) | Physical constants and altitude-dependent parameter arrays. |
| [main](main.md) | Command-line entry point tying everything together. |
| [grid_cache](grid_cache.md) | Solves each grid once and caches it as CSV in `data/`. |
| [methods](methods.md) | Plotting methods (growth, frequency, ratio, phase velocity). |
| [csvmethods](csvmethods.md) | Post-processing of cached CSVs (resonance bands, a_i0 differences). |
| [replot_variants](replot_variants.md) | Batch script regenerating all variant grids/plots (normal, a_i0 = 0, doubled E-field). |
| [compare_electrons](compare_electrons.md) | Script comparing fluid vs. fully kinetic electrons on the same k-grid. |