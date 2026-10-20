# checkpoint_dispersion.py
#
# Run:
#     python checkpoint_dispersion.py
#
# This script independently validates candidates. A solver's success flag
# or a small absolute residual is NOT sufficient to accept a root.

import inspect
import warnings

import numpy as np
from scipy.optimize import least_squares

import constants as c
import eq70 as model
import root as root_module


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

RUN = 5                          # 104 km
TEST_POINTS = [
    (1.0, 0.0),
    (1.0, -0.5),
]

REL_TOL = 1e-8
DEDUP_TOL = 1e-5                 # Relative to the frequency scale
BRANCH_MARGIN = 1e-6             # Relative to nu_i

# The search windows are numerical diagnostics, not physical validity limits.
# Enlarge them to check whether additional roots appear.
WINDOWS = [
    (-4.0, 4.0, 2.0),            # xmin, xmax, ymax in scaled frequency
    (-8.0, 8.0, 4.0),
]

CSV_FILE = "checkpoint_roots.csv"


# ------------------------------------------------------------
# Model adapters
# ------------------------------------------------------------

def supports_prefactor(func):
    return "ion_prefactor" in inspect.signature(func).parameters


TERMS_HAS_PREFACTOR = supports_prefactor(model.dispersion_terms)
WRAPPER_HAS_PREFACTOR = supports_prefactor(model.eq_70)

MODES = ["A83", "printed73"] if TERMS_HAS_PREFACTOR else [None]


def terms(omega, kx, ky, mode=None):
    kwargs = {"run": RUN}
    if mode is not None and TERMS_HAS_PREFACTOR:
        kwargs["ion_prefactor"] = mode

    result = model.dispersion_terms(omega, kx, ky, **kwargs)
    if len(result) != 3:
        raise ValueError("dispersion_terms must return (I_i, I_e, Q).")

    return tuple(complex(value) for value in result)


def D(omega, kx, ky, mode=None):
    Ii, Ie, Q = terms(omega, kx, ky, mode)
    return Ii + Ie + Q


def wrapper_D(omega, kx, ky, mode=None):
    kwargs = {"run": RUN}
    if mode is not None and WRAPPER_HAS_PREFACTOR:
        kwargs["ion_prefactor"] = mode
    return model.eq_70(omega, kx, ky, **kwargs)


def relative_residual(omega, kx, ky, mode=None):
    try:
        Ii, Ie, Q = terms(omega, kx, ky, mode)
    except (FloatingPointError, OverflowError, ZeroDivisionError):
        return np.inf

    if not np.all(np.isfinite([Ii, Ie, Q])):
        return np.inf

    scale = abs(Ii) + abs(Ie) + abs(Q)
    if scale == 0:
        return np.inf

    return abs(Ii + Ie + Q) / scale


def frequency_scale(kx, ky):
    return max(float(c.nu_i[RUN]), np.hypot(kx, ky) * abs(c.V_d))


def in_upper_branch(omega):
    nui = float(c.nu_i[RUN])
    return omega.imag > -nui + BRANCH_MARGIN * nui


def accepted(omega, kx, ky, mode=None):
    return (
        np.isfinite(omega)
        and in_upper_branch(omega)
        and relative_residual(omega, kx, ky, mode) < REL_TOL
    )


def deduplicate(candidates, scale):
    unique = []

    for candidate in candidates:
        if not any(
            abs(candidate - existing) < DEDUP_TOL * scale
            for existing in unique
        ):
            unique.append(candidate)

    return sorted(unique, key=lambda value: value.imag, reverse=True)


# ------------------------------------------------------------
# Step 2: wrapper versus term sum
# ------------------------------------------------------------

def check_wrapper():
    print("\n=== STEP 2: wrapper versus term sum ===")

    # If the wrapper has no mode argument, compare its default configuration.
    modes = MODES if WRAPPER_HAS_PREFACTOR else [None]

    for mode in modes:
        for kx, ky in TEST_POINTS:
            scale = frequency_scale(kx, ky)

            for omega in [
                0.2 * scale + 0.05j * scale,
                -0.4 * scale + 0.1j * scale,
                0.5 * scale - 0.25j * c.nu_i[RUN],
            ]:
                summed = D(omega, kx, ky, mode)
                wrapped = wrapper_D(omega, kx, ky, mode)

                assert np.isfinite(summed), (
                    f"Nonfinite term sum at omega={omega}, k=({kx}, {ky})"
                )
                assert np.isfinite(wrapped), (
                    f"Nonfinite wrapper at omega={omega}, k=({kx}, {ky})"
                )

                np.testing.assert_allclose(
                    wrapped, summed, rtol=1e-12, atol=1e-14
                )

        print(f"PASS: mode={mode or 'default'}")

    print(
        "Note: this checks wrapper consistency, not whether the "
        "refactor preserved your old formula."
    )


# ------------------------------------------------------------
# Step 3: distant frequencies must not pass root validation
# ------------------------------------------------------------

def check_remote_candidates():
    print("\n=== STEP 3: remote candidate rejection ===")

    for mode in MODES:
        for kx, ky in TEST_POINTS:
            scale = frequency_scale(kx, ky)

            # Near-imaginary directions can produce very different Z behavior.
            # This direction is chosen as a diagnostic, not a universal
            # asymptotic test for every possible model parameter.
            for factor in [1e2, 1e4]:
                omega = factor * scale * (1.0 + 0.1j)

                residual = relative_residual(omega, kx, ky, mode)

                with np.errstate(all="ignore"):
                    absolute = abs(D(omega, kx, ky, mode))

                print(
                    f"mode={mode or 'default':>9} "
                    f"k=({kx: .2f},{ky: .2f}) "
                    f"|omega|={abs(omega):.3e} "
                    f"|D|={absolute:.3e} "
                    f"relative={residual:.3e}"
                )

                assert not accepted(omega, kx, ky, mode), (
                    "A deliberately remote candidate passed validation. "
                    "Inspect the dispersion terms before continuing."
                )

    print("PASS: remote test candidates were rejected.")


# ------------------------------------------------------------
# Step 4: inspect existing Newton solver
# ------------------------------------------------------------

def extract_newton_result(result):
    # Supports your original complex return and the proposed dictionary return.
    if isinstance(result, dict):
        return (
            complex(result["omega"]),
            result.get("converged", "not reported"),
            result.get("reason", "not reported"),
        )

    return complex(result), "not reported", "legacy return"


def check_newton():
    print("\n=== STEP 4: Newton candidate validation ===")

    solver = getattr(
        root_module, "custom_newton_raphson_2d", None
    )

    if solver is None:
        print("SKIP: custom_newton_raphson_2d was not found.")
        return

    params = inspect.signature(solver).parameters

    for mode in MODES:
        for kx, ky in TEST_POINTS:
            scale = frequency_scale(kx, ky)
            seeds = list(
                root_module.fluid_seed_pair_2d(kx, ky, run=RUN)
            )
            seeds.append(0.5 * scale + 0.1j * scale)

            for seed in seeds:
                if not in_upper_branch(seed):
                    print(
                        f"SKIP seed outside upper branch: {seed:.6g}"
                    )
                    continue

                func = lambda w, x, y: D(w, x, y, mode)
                validator = lambda w, x, y: relative_residual(
                    w, x, y, mode
                )

                kwargs = {
                    "func": func,
                    "initial_guess": seed,
                    "kx": kx,
                    "ky": ky,
                }

                # Recognize either of these callback names if you added one.
                if "relative_residual_func" in params:
                    kwargs["relative_residual_func"] = validator
                elif "residual_func" in params:
                    kwargs["residual_func"] = validator

                try:
                    result = solver(**kwargs)
                    omega, reported, reason = extract_newton_result(result)
                except Exception as exc:
                    print(f"Newton raised {type(exc).__name__}: {exc}")
                    continue

                residual = relative_residual(omega, kx, ky, mode)
                valid = accepted(omega, kx, ky, mode)

                print(
                    f"mode={mode or 'default'} "
                    f"k=({kx},{ky}) seed={seed:.6g}\n"
                    f"  candidate={omega:.9g} "
                    f"relative={residual:.3e} "
                    f"validated={valid}\n"
                    f"  solver_converged={reported}; reason={reason}"
                )

                if reported is True and not valid:
                    warnings.warn(
                        "Newton reported convergence, but independent "
                        "validation rejected the candidate."
                    )

    print(
        "Newton failures here are diagnostic results, not automatically "
        "script failures."
    )


# ------------------------------------------------------------
# Step 5: mean ion-temperature bookkeeping
# ------------------------------------------------------------

def check_temperature():
    print("\n=== STEP 5: ion-temperature balance ===")

    ai = c.e * c.E_0 / c.m_i
    Ui = ai / np.asarray(c.nu_i)

    expected = c.T_n + c.m_i * Ui**2 / 3.0

    # Your coordinate convention: y is parallel to E0;
    # x and the omitted z axis are equivalent transverse axes.
    calculated = c.m_i * (
        2.0 * np.asarray(c.V_Tx)**2
        + np.asarray(c.V_Ty)**2
    ) / 3.0

    np.testing.assert_allclose(
        calculated, expected, rtol=1e-12, atol=0
    )
    print("PASS: thermal velocities satisfy the mean energy balance.")

    if hasattr(c, "T_eff_i"):
        stored = np.asarray(c.T_eff_i)
        error = np.max(np.abs(stored - expected) / expected)

        if error < 1e-12:
            print("PASS: stored T_eff_i is consistent.")
        else:
            warnings.warn(
                f"Stored T_eff_i has maximum relative error {error:.3e}. "
                "Use m_i*(2*V_Tx**2 + V_Ty**2)/3."
            )

    if not TERMS_HAS_PREFACTOR:
        print(
            "SKIP prefactor comparison: add an explicit ion_prefactor "
            "argument to dispersion_terms."
        )
    else:
        print("Both prefactor modes will be searched below.")

    print(
        "The isotropic/no-field equality test requires explicit ion-parameter "
        "overrides in your model. It is not tested by this script."
    )


# ------------------------------------------------------------
# Steps 6–7: bounded, scaled single-point searches
# ------------------------------------------------------------

def bounded_search(kx, ky, mode, window):
    scale = frequency_scale(kx, ky)
    nui = float(c.nu_i[RUN])

    xmin, xmax, ymax = window
    ymin = (-nui + BRANCH_MARGIN * nui) / scale

    lower = np.array([xmin, ymin])
    upper = np.array([xmax, ymax])

    # Constant residual scaling. Do not normalize by |Ii|+|Ie|+|Q|
    # inside the solver: use that only for independent validation.
    reference = terms(0.2 * scale + 0.1j * scale, kx, ky, mode)
    residual_scale = max(sum(abs(value) for value in reference), 1e-12)

    def objective(x):
        omega = scale * complex(x[0], x[1])

        with np.errstate(all="ignore"):
            value = D(omega, kx, ky, mode) / residual_scale

        if not np.isfinite(value):
            # Finite penalty for an unusable evaluation.
            # Such candidates are still rejected independently.
            return np.array([1e12, 1e12])

        return np.array([value.real, value.imag])

    # Deterministic seed grid, plus the corrected fluid seeds.
    xs = np.linspace(xmin + 0.05, xmax - 0.05, 7)
    ys = np.array([
        0.75 * ymin,
        0.0,
        0.25 * ymax,
        0.75 * ymax,
    ])

    seeds = [
        scale * complex(x, y)
        for x in xs
        for y in ys
    ]
    seeds.extend(
        root_module.fluid_seed_pair_2d(kx, ky, run=RUN)
    )

    candidates = []
    attempts = 0

    for seed in seeds:
        if not in_upper_branch(seed):
            continue

        x0 = np.array([seed.real / scale, seed.imag / scale])
        if np.any(x0 <= lower) or np.any(x0 >= upper):
            continue

        attempts += 1

        try:
            result = least_squares(
                objective,
                x0,
                bounds=(lower, upper),
                xtol=1e-11,
                ftol=1e-11,
                gtol=1e-11,
                max_nfev=500,
            )
        except (ValueError, FloatingPointError, OverflowError):
            continue

        omega = scale * complex(*result.x)

        # Acceptance is independent of result.success.
        if accepted(omega, kx, ky, mode):
            candidates.append(omega)

    return deduplicate(candidates, scale), attempts


def check_single_point_roots():
    print("\n=== STEPS 6–7: bounded single-point searches ===")

    csv_rows = []

    for mode in MODES:
        for kx, ky in TEST_POINTS:
            previous = []

            for window_index, window in enumerate(WINDOWS):
                found, attempts = bounded_search(kx, ky, mode, window)

                print(
                    f"\nmode={mode or 'default'}, k=({kx},{ky}), "
                    f"window={window}, attempts={attempts}, "
                    f"validated roots={len(found)}"
                )

                for omega in found:
                    Ii, Ie, Q = terms(omega, kx, ky, mode)
                    residual = relative_residual(omega, kx, ky, mode)

                    print(
                        f"  omega={omega.real:+.9e}"
                        f"{omega.imag:+.9e}j  relative={residual:.3e}\n"
                        f"    Ii={Ii:.6g}, Ie={Ie:.6g}, Q={Q:.6g}"
                    )

                    csv_rows.append([
                        mode or "default",
                        window_index,
                        kx,
                        ky,
                        omega.real,
                        omega.imag,
                        residual,
                    ])

                scale = frequency_scale(kx, ky)
                lost = [
                    old for old in previous
                    if not any(
                        abs(old - new) < DEDUP_TOL * scale
                        for new in found
                    )
                ]

                if lost:
                    warnings.warn(
                        "The larger window did not recover some earlier roots. "
                        "Increase seed coverage; search robustness is unresolved."
                    )

                previous = found

            if not previous:
                warnings.warn(
                    f"No validated roots found for mode={mode}, "
                    f"k=({kx},{ky}). This does not prove stability."
                )

    with open(CSV_FILE, "w", encoding="utf-8") as file:
        file.write("mode,window,kx,ky,omega_real,omega_imag,relative_residual\n")
        for row in csv_rows:
            file.write(",".join(map(str, row)) + "\n")

    print(f"\nSaved validated candidates to {CSV_FILE}")
    print("These searches do not establish root completeness.")


# ------------------------------------------------------------
# Step 8: missing-data convention
# ------------------------------------------------------------

def check_missing_data():
    print("\n=== STEP 8: missing-data convention ===")

    roots = np.full((2, 2), complex(np.nan, np.nan))
    growth = np.full((2, 2), np.nan)

    roots[0, 0] = 100.0 + 20.0j
    growth[0, 0] = roots[0, 0].imag

    neighbor_seeds = [
        value for value in roots.ravel()
        if np.isfinite(value)
    ]

    assert len(neighbor_seeds) == 1
    assert np.isnan(growth[1, 1])
    assert not np.isfinite(roots[1, 1])

    print("PASS: demonstration preserves failed points as NaN.")
    print(
        "This tests the convention only. Update contour_plot to use the "
        "same initialization and finite-neighbor checks."
    )


def main():
    print(f"Checkpoint run={RUN}, altitude={c.alt_eff[RUN]} km")
    print(f"nu_i={c.nu_i[RUN]:.6g} s^-1")
    print(f"Relative acceptance tolerance={REL_TOL:.1e}")

    check_wrapper()
    check_remote_candidates()
    check_newton()
    check_temperature()
    check_single_point_roots()
    check_missing_data()

    print("\nCheckpoint diagnostics completed.")


if __name__ == "__main__":
    main()