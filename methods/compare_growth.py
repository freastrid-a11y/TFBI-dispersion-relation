import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_growth_csv(path, growth_name):
    """Read a growth-rate CSV and keep the kx, ky, and growth-rate columns."""
    df = pd.read_csv(path)

    required = {"kx_m^-1", "ky_m^-1", "growth_rate_s^-1"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{path} is missing columns: {', '.join(sorted(missing))}"
        )

    df = df[list(required)].copy()
    df = df.rename(columns={"growth_rate_s^-1": growth_name})
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Compare growth rates from two CSV files and plot their relative difference."
    )
    parser.add_argument("file1", help="CSV containing growth rate g1")
    parser.add_argument("file2", help="CSV containing growth rate g2")
    parser.add_argument(
        "--output",
        default="growth_rate_comparison.csv",
        help="Output CSV filename (default: growth_rate_comparison.csv)",
    )
    args = parser.parse_args()

    df1 = read_growth_csv(args.file1, "g1")
    df2 = read_growth_csv(args.file2, "g2")

    # Match points with the same kx and ky values.
    merged = df1.merge(
        df2,
        on=["kx_m^-1", "ky_m^-1"],
        how="inner",
        validate="one_to_one",
    )

    if merged.empty:
        raise ValueError(
            "The files have no matching (kx, ky) points. "
            "Check that both CSVs use the same wave-vector grid."
        )

    if len(merged) < len(df1) or len(merged) < len(df2):
        print(
            f"Warning: only {len(merged)} matching (kx, ky) points; "
            "some rows were omitted."
        )

    denominator = (merged["g1"] + merged["g2"]) / 2.0

    # The expression is undefined when the average growth rate is zero.
    merged["relative_difference"] = np.where(
        denominator != 0.0,
        (merged["g1"] - merged["g2"]) / denominator,
        np.nan,
    )

    merged.to_csv(args.output, index=False)
    print(f"Saved comparison data to {args.output}")

    # Arrange the values on a 2-D grid for plotting.
    grid = merged.pivot(
        index="ky_m^-1",
        columns="kx_m^-1",
        values="relative_difference",
    ).sort_index().sort_index(axis=1)

    kx = grid.columns.to_numpy(dtype=float)
    ky = grid.index.to_numpy(dtype=float)
    values = grid.to_numpy(dtype=float)

    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError("No finite relative-difference values to plot.")

    # Use a color scale centered at zero.
    limit = np.max(np.abs(finite))
    if limit == 0:
        limit = 1.0

    fig, ax = plt.subplots(figsize=(8, 6))
    mesh = ax.pcolormesh(
        kx,
        ky,
        values,
        shading="auto",
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
    )

    fig.colorbar(
        mesh,
        ax=ax,
        label=r"$(g_1-g_2)/[(g_1+g_2)/2]$",
    )
    ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
    ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    ax.set_title("Relative difference in growth rate")
    fig.tight_layout()
    plt.savefig("Growth_rate_difference_110km.png")
    plt.show()


if __name__ == "__main__":
    main()