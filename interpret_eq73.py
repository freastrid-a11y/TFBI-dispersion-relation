import numpy as np
import matplotlib.pyplot as plt

from eq73 import eq_73
from root import find_roots
import constants as C


# ============================================================
# USER SETTINGS
# ============================================================

RUN = 5                  # alt_eff[5] = 104 km

K_MIN = 0.05             # Nonzero: avoid the singular k = 0 point
N_SMALL = 14             # Logarithmic samples for small-k diagnostics
N_LARGE = 50             # Samples through the plotted wave-number range

# Limits of your original 2D plot
KX_LIMIT = 15.0
KY_LIMIT = 25.0

# If you know the correct initial root, supply it here, in rad/s.
# Keys are the direction names below.
#
# Example:
# INITIAL_ROOTS = {"x direction": 10.0 + 0.1j}
INITIAL_ROOTS = {}

# A large predictor error indicates that branch tracking needs inspection.
JUMP_WARNING = 0.35


# ============================================================
# EQUILIBRIUM DRIFTS
# ============================================================

def equilibrium_drifts(run):
    """Use the same equilibrium velocities as eq73.py."""
    ke = C.kappa_e[run]

    electron = np.array([
        C.V_d * ke**2 / (1.0 + ke**2),
        -C.V_d * ke / (1.0 + ke**2),
    ])

    ion = np.array([
        0.0,
        C.e * C.E_0 / (C.m_i * C.nu_i[run]),
    ])

    return electron, ion


def approximate_drift_frequency(k_vector, run):
    """
    Approximate low-frequency drift frequency used ONLY to choose
    the first root. It is not an exact solution of Eq. 73.
    """
    ve, ui = equilibrium_drifts(run)
    psi = C.Psi_0[run]

    return float(
        np.dot(k_vector, ve + psi * ui) / (1.0 + psi)
    )


def maximum_k_in_rectangle(direction):
    """Stop each cut at the boundary of the original plotting domain."""
    candidates = []

    if abs(direction[0]) > 1e-12:
        candidates.append(KX_LIMIT / abs(direction[0]))

    if abs(direction[1]) > 1e-12:
        candidates.append(KY_LIMIT / abs(direction[1]))

    return min(candidates)


def make_k_values(k_max):
    small_end = min(0.8, 0.25 * k_max)

    small = np.geomspace(K_MIN, small_end, N_SMALL)
    large = np.linspace(small_end, k_max, N_LARGE)

    return np.unique(np.concatenate([small, large]))


# ============================================================
# ROOT SELECTION AND CONTINUATION
# ============================================================

def candidate_roots(k_vector, run, seeds):
    """Multistart search plus continuation seeds."""
    return find_roots(
        eq_73,
        float(k_vector[0]),
        float(k_vector[1]),
        run,
        guesses=seeds,
        include_default_guesses=True,
    )


def continue_cut(direction, k_values, run, initial_root=None):
    """
    Follow a root along k_vector = k * direction.

    Selection:
      - first point: closest root to an approximate drift frequency
      - later points: closest root to a continuation predictor

    The independently highest-growth candidate is retained for comparison.
    """
    n = len(k_values)

    tracked = np.full(n, np.nan + 1j * np.nan, dtype=complex)
    highest = np.full(n, np.nan + 1j * np.nan, dtype=complex)

    raw_residual = np.full(n, np.nan)
    predictor_error = np.full(n, np.nan)
    candidate_count = np.zeros(n, dtype=int)

    history = []

    for j, k in enumerate(k_values):
        kv = k * direction

        if not history:
            target = (
                complex(initial_root)
                if initial_root is not None
                else complex(approximate_drift_frequency(kv, run))
            )

        elif len(history) == 1:
            # Expected leading-order drift scaling: omega ~ k.
            last = history[-1]
            target = tracked[last] * k / k_values[last]

        else:
            # Linear extrapolation from the last two accepted points.
            previous, last = history[-2:]

            slope = (
                (tracked[last] - tracked[previous])
                / (k_values[last] - k_values[previous])
            )

            target = tracked[last] + slope * (k - k_values[last])

        seeds = [target]

        if history:
            seeds.append(tracked[history[-1]])

        roots = candidate_roots(kv, run, seeds)
        candidate_count[j] = len(roots)

        if not roots:
            # Leave a gap rather than inventing a root.
            continue

        highest[j] = max(roots, key=lambda w: w.imag)
        selected = min(roots, key=lambda w: abs(w - target))

        tracked[j] = selected

        comparison_scale = max(
            1.0,
            abs(target),
            abs(selected),
            np.linalg.norm(kv)
            * np.sqrt(C.C_s_squared),
        )

        predictor_error[j] = (
            abs(selected - target) / comparison_scale
        )

        raw_residual[j] = abs(
            eq_73(selected, kv[0], kv[1], run)
        )

        history.append(j)

    return {
        "k": k_values,
        "omega": tracked,
        "highest": highest,
        "residual": raw_residual,
        "predictor_error": predictor_error,
        "candidate_count": candidate_count,
    }


# ============================================================
# INDEPENDENT REVERSAL-SYMMETRY CHECK
# ============================================================

def reversal_check(direction, data, run):
    """
    Independently solve at -k, then compare with -conjugate(omega(k)).

    An unseeded multistart search is used so that agreement is not
    imposed by initializing the solver at the symmetry-predicted root.
    """
    valid = np.flatnonzero(np.isfinite(data["omega"]))

    if len(valid) == 0:
        return []

    sample_indices = np.unique(
        valid[np.linspace(0, len(valid) - 1, 6).astype(int)]
    )

    results = []

    for j in sample_indices:
        k = data["k"][j]
        omega = data["omega"][j]
        kv_negative = -k * direction

        roots = candidate_roots(kv_negative, run, seeds=[])

        if not roots:
            results.append((k, np.nan))
            continue

        expected = -np.conjugate(omega)
        matched = min(roots, key=lambda w: abs(w - expected))

        relative_error = (
            abs(matched - expected) / max(1.0, abs(expected))
        )

        results.append((k, relative_error))

    return results


# ============================================================
# PRINTED DIAGNOSTICS
# ============================================================

def print_summary(name, direction, data, run):
    k = data["k"]
    omega = data["omega"]
    highest = data["highest"]

    valid = np.isfinite(omega)
    indices = np.flatnonzero(valid)

    print("\n" + "=" * 72)
    print(name)
    print(f"Direction = ({direction[0]:+.5f}, {direction[1]:+.5f})")
    print(f"Accepted tracked roots: {len(indices)}/{len(k)}")

    if len(indices) == 0:
        print("No branch could be tracked.")
        return

    print("\nFirst small-k points:")
    print(" k [1/m]      Re(omega)       Im(omega)       Re/k       Im/k^2")

    for j in indices[:6]:
        print(
            f"{k[j]:9.4g}  "
            f"{omega[j].real:13.5g}  "
            f"{omega[j].imag:13.5g}  "
            f"{omega[j].real / k[j]:10.5g}  "
            f"{omega[j].imag / k[j]**2:10.5g}"
        )

    print(
        "\nLargest raw |eq_73| among tracked points: "
        f"{np.nanmax(data['residual']):.3e}"
    )
    print(
        "Note: this raw residual belongs to your scaled dispersion equation;"
        "\nits magnitude alone is not a physical accuracy estimate."
    )

    j_peak = indices[np.argmax(omega[indices].imag)]

    print(
        f"\nLargest tracked growth = {omega[j_peak].imag:.6g} s^-1"
        f" at k = {k[j_peak]:.6g} m^-1"
    )

    if j_peak == len(k) - 1:
        print("WARNING: the tracked growth maximum is at the domain boundary.")

    comparison_valid = valid & np.isfinite(highest)

    n_different = np.count_nonzero(
        comparison_valid
        & (
            np.abs(highest - omega)
            > 1e-4 * np.maximum(1.0, np.abs(omega))
        )
    )

    print(
        f"Highest-growth selection differs from tracked branch at "
        f"{n_different} points."
    )

    # The first point is selected, rather than continuation-predicted.
    jump_indices = np.flatnonzero(
        (data["predictor_error"] > JUMP_WARNING)
        & (np.arange(len(k)) > 0)
    )

    if len(jump_indices):
        print(
            "WARNING: large continuation predictor errors at k =",
            np.round(k[jump_indices], 4),
        )
        print("Refine sampling and inspect possible branch switches.")

    print("\nIndependent wave-vector reversal check:")
    for kval, error in reversal_check(direction, data, run):
        if np.isfinite(error):
            print(f"  k = {kval:8.4g}: relative mismatch = {error:.3e}")
        else:
            print(f"  k = {kval:8.4g}: no negative-k root found")


# ============================================================
# PLOTTING
# ============================================================

def plot_branch_diagnostics(results, run):
    fig, axes = plt.subplots(
        len(results),
        4,
        figsize=(18, 4 * len(results)),
        squeeze=False,
    )

    for row, (name, direction, data) in enumerate(results):
        k = data["k"]
        omega = data["omega"]
        highest = data["highest"]

        ax = axes[row, 0]
        ax.plot(k, omega.real, label="Tracked branch")
        ax.plot(
            k,
            highest.real,
            "--",
            alpha=0.75,
            label="Highest-growth candidate",
        )
        ax.set_ylabel(f"{name}\nRe(ω) [rad/s]")
        ax.set_title("Real frequency")
        ax.legend(fontsize=8)

        ax = axes[row, 1]
        ax.plot(k, omega.imag, label="Tracked branch")
        ax.plot(
            k,
            highest.imag,
            "--",
            alpha=0.75,
            label="Highest growth found",
        )
        ax.axhline(0, color="black", lw=0.8)

        suspicious = (
            data["predictor_error"] > JUMP_WARNING
        ) & (np.arange(len(k)) > 0)

        ax.scatter(
            k[suspicious],
            omega.imag[suspicious],
            color="red",
            marker="x",
            label="Large predictor error",
            zorder=5,
        )

        ax.set_ylabel("Im(ω) [s⁻¹]")
        ax.set_title("Growth / damping")
        ax.legend(fontsize=8)

        # Show only the small-k portion for scaling checks.
        small = k <= min(0.8, k[-1])

        ax = axes[row, 2]
        ax.semilogx(
            k[small],
            omega.real[small] / k[small],
            "o-",
        )
        ax.set_ylabel("Re(ω) / k [m/s]")
        ax.set_title("Small-k phase-speed check")

        ax = axes[row, 3]
        ax.semilogx(
            k[small],
            omega.imag[small] / k[small]**2,
            "o-",
        )
        ax.axhline(0, color="black", lw=0.8)
        ax.set_ylabel("Im(ω) / k² [m²/s]")
        ax.set_title("Small-k growth-scaling check")

        for ax in axes[row]:
            ax.set_xlabel("|k| [m⁻¹]")
            ax.grid(alpha=0.25)

    fig.suptitle(
        f"Eq. 73 branch diagnostics — {C.alt_eff[run]:g} km",
        fontsize=15,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])


def plot_frequency_ratios(results, run):
    """
    These ratios are diagnostics, not a complete validity test.
    Their relevance depends on the assumptions behind Eq. 73.
    """
    ve, _ = equilibrium_drifts(run)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    for name, direction, data in results:
        k = data["k"]
        omega = data["omega"]

        ion_ratio = np.abs(omega) / C.nu_i[run]

        electron_frame_frequency = (
            omega - k * np.dot(direction, ve)
        )

        electron_ratio = (
            np.abs(electron_frame_frequency) / C.nu_e[run]
        )

        axes[0].loglog(k, ion_ratio, label=name)
        axes[1].loglog(k, electron_ratio, label=name)

    axes[0].set_title("Ion-frame frequency / ion collisions")
    axes[0].set_ylabel("|ω| / νᵢ")

    axes[1].set_title("Electron-frame frequency / electron collisions")
    axes[1].set_ylabel("|ω − k·Vₑ| / νₑ")

    for ax in axes:
        ax.axhline(1.0, color="black", ls="--", label="Ratio = 1")
        ax.set_xlabel("|k| [m⁻¹]")
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=8)

    fig.suptitle(
        "Frequency ordering diagnostics — not automatic pass/fail tests"
    )
    fig.tight_layout()


# ============================================================
# MAIN
# ============================================================

def main():
    ve, ui = equilibrium_drifts(RUN)
    relative = ve - ui
    relative_direction = relative / np.linalg.norm(relative)

    directions = [
        ("x direction", np.array([1.0, 0.0])),
        ("y direction", np.array([0.0, 1.0])),
        ("relative drift", relative_direction),
    ]

    print(f"Altitude: {C.alt_eff[RUN]:g} km")
    print(f"Electron drift: {ve} m/s")
    print(f"Ion drift:      {ui} m/s")
    print(f"Relative drift: {relative} m/s")
    print(f"Ion collisions:      {C.nu_i[RUN]:.6g} s^-1")
    print(f"Electron collisions: {C.nu_e[RUN]:.6g} s^-1")

    results = []

    for name, direction in directions:
        k_max = maximum_k_in_rectangle(direction)
        k_values = make_k_values(k_max)

        print(f"\nTracking {name} ...")

        data = continue_cut(
            direction,
            k_values,
            RUN,
            initial_root=INITIAL_ROOTS.get(name),
        )

        print_summary(name, direction, data, RUN)
        results.append((name, direction, data))

    plot_branch_diagnostics(results, RUN)
    plt.savefig("diagnostics")
    plot_frequency_ratios(results, RUN)
    plt.savefig("frequency_ratios")
    plt.show()


if __name__ == "__main__":
    main()