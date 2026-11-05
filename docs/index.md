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
| [eq73](eq73.md) | The TFBI dispersion relation (manuscript Eq. 73) — the core math. |
| [root](root.md) | Modified plasma dispersion function, multistart root finder, grid sweep. |
| [constants](constants.md) | Physical constants and altitude-dependent parameter arrays. |
| [main](main.md) | Command-line entry point tying everything together. |
| [grid_cache](grid_cache.md) | Solves each grid once and caches it as CSV in `data/`. |
| [methods](methods.md) | Plotting methods (growth, frequency, ratio, phase velocity). |
| [csvmethods](csvmethods.md) | Post-processing of cached CSVs (resonance bands, a_i0 differences). |