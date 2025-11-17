#!/usr/bin/env python3
"""
Analyze weight sweep results and visualize correlations
"""

import pandas as pd
import sys


def main():
    # Read results
    try:
        df = pd.read_csv("weight_sweep_results.csv")
    except FileNotFoundError:
        print("ERROR: weight_sweep_results.csv not found!")
        sys.exit(1)

    print("=" * 70)
    print(" WEIGHT SWEEP ANALYSIS")
    print("=" * 70)
    print()

    # Summary statistics
    print(f"Total configurations tested: {len(df)}")
    print(f"Best DockQ achieved: {df['best_dockq'].max():.3f}")
    print(f"Worst DockQ: {df['best_dockq'].min():.3f}")
    print(f"Average DockQ: {df['best_dockq'].mean():.3f}")
    print()

    # Top 10 configurations
    print("=" * 70)
    print("TOP 10 CONFIGURATIONS BY DOCKQ")
    print("=" * 70)
    print()

    df_sorted = df.sort_values("best_dockq", ascending=False)

    print(
        f"{'Rank':<5} {'DockQ':<8} {'VDW':<6} {'Elec':<6} {'Desolv':<8} {'AIR':<8} {'RMSD':<7} {'Conv':<5}"
    )
    print("-" * 70)

    for i, (idx, row) in enumerate(df_sorted.head(10).iterrows(), 1):
        print(
            f"{i:<5} {row['best_dockq']:<8.3f} {row['w_vdw']:<6.1f} {row['w_elec']:<6.1f} "
            f"{row['w_desolv']:<8.1f} {row['w_air']:<8.1f} {row['best_rmsd']:<7.2f} {row['convergence_gen']:<5}"
        )

    print()

    # Correlation analysis
    print("=" * 70)
    print(" WEIGHT IMPACT ANALYSIS")
    print("=" * 70)
    print()

    print("Correlation with DockQ:")
    for weight in ["w_vdw", "w_elec", "w_desolv", "w_air"]:
        corr = df[weight].corr(df["best_dockq"])
        direction = "↑" if corr > 0 else "↓"
        strength = (
            "Strong" if abs(corr) > 0.5 else "Moderate" if abs(corr) > 0.3 else "Weak"
        )
        print(f"  {weight:12} {direction} {corr:+.3f} ({strength})")

    print()

    # Key insights
    print("=" * 70)
    print(" RESULTS")
    print("=" * 70)
    print()

    # Best configuration
    best_row = df_sorted.iloc[0]
    print("✓ Best configuration:")
    print(
        f"  VDW={best_row['w_vdw']:.1f}, Elec={best_row['w_elec']:.1f}, "
        f"Desolv={best_row['w_desolv']:.1f}, AIR={best_row['w_air']:.1f}"
    )
    print(f"  → DockQ={best_row['best_dockq']:.3f}")
    print()

    # Baseline comparison
    # HACK: Baseline is hardcoded
    baseline = df[
        (df["w_vdw"] == 1.0)
        & (df["w_elec"] == 0.5)
        & (df["w_desolv"] == 0.5)
        & (df["w_air"] == 100.0)
    ]

    if not baseline.empty:
        baseline_dockq = baseline["best_dockq"].values[0]
        improvement = ((best_row["best_dockq"] - baseline_dockq) / baseline_dockq) * 100
        print("✓ Improvement over baseline:")
        print(f"  Baseline DockQ: {baseline_dockq:.3f}")
        print(f"  Best DockQ: {best_row['best_dockq']:.3f}")
        print(f"  Improvement: {improvement:+.1f}%")
        print()

    # AIR weight analysis
    print("✓ AIR weight effect:")
    air_groups = df.groupby("w_air")["best_dockq"].agg(["mean", "max", "count"])
    for air_val, stats in air_groups.iterrows():
        print(
            f"  AIR={air_val:6.1f}: avg DockQ={stats['mean']:.3f}, "
            f"max={stats['max']:.3f} (n={int(stats['count'])})"
        )
    print()


if __name__ == "__main__":
    main()
