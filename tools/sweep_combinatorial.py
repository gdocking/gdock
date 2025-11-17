#!/usr/bin/env python3
"""
Combinatorial grid search for weight optimization
"""

# Yeah, python (:

import subprocess
import pandas as pd
import re
import sys
from pathlib import Path
from datetime import datetime
import itertools

# Define grid values for each parameter
GRID = {
    "w_vdw": [0.1, 0.5, 1.0, 2.0, 3.0],
    "w_elec": [0.1, 0.5, 1.0, 1.5],
    "w_desolv": [0.1, 0.5, 1.0, 1.5],
    "w_air": [1.0, 10.0, 50.0, 100.0, 150.0, 200.0, 500.0, 1000.0],
}


def build():
    """Build gdock in release mode"""
    subprocess.run(
        ["cargo", "build", "--release"], capture_output=True, text=True, check=True
    )


def run(w_vdw, w_elec, w_desolv, w_air, log_file):
    """Run gdock with specified weights"""
    cmd = [
        "./target/release/gdock",
        "run",
        "--w_vdw",
        str(w_vdw),
        "--w_elec",
        str(w_elec),
        "--w_desolv",
        str(w_desolv),
        "--w_air",
        str(w_air),
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )

        # Write output to log file
        log_file.write_text(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print("  ERROR: Run failed!")
        print(e.stderr)
        return False


def extract_best_dockq(log_file):
    """Extract the best DockQ value from all generations"""
    try:
        content = log_file.read_text()

        # Find all DockQ values in generation lines
        dockq_pattern = r"DockQ=(\d+\.\d+)"
        matches = re.findall(dockq_pattern, content)

        if matches:
            dockq_values = [float(x) for x in matches]
            return max(dockq_values)
        else:
            return 0.0
    except Exception as e:
        print(f"  WARNING: Could not extract DockQ: {e}")
        return 0.0


def extract_final_metrics(log_file):
    """Extract metrics from the final generation"""
    try:
        content = log_file.read_text()

        # Get the last generation line
        gen_lines = [line for line in content.split("\n") if line.startswith("Gen #")]
        if not gen_lines:
            return None, None, None

        last_line = gen_lines[-1]

        # Extract score, RMSD, FNAT
        score_match = re.search(r"Score=(\d+\.\d+)", last_line)
        rmsd_match = re.search(r"RMSD=(\d+\.\d+)", last_line)
        fnat_match = re.search(r"FNAT=(\d+\.\d+)", last_line)

        score = float(score_match.group(1)) if score_match else None
        rmsd = float(rmsd_match.group(1)) if rmsd_match else None
        fnat = float(fnat_match.group(1)) if fnat_match else None

        return score, rmsd, fnat
    except Exception as e:
        print(f"  WARNING: Could not extract final metrics: {e}")
        return None, None, None


def extract_convergence_gen(log_file):
    """Extract convergence generation"""
    try:
        content = log_file.read_text()

        match = re.search(r"Converged at generation (\d+)", content)
        if match:
            return int(match.group(1))
        else:
            return None
    except Exception as _:
        return None


def generate_configurations():
    """Generate all combinations of weight values"""
    # Create all combinations
    combinations = list(
        itertools.product(
            GRID["w_vdw"], GRID["w_elec"], GRID["w_desolv"], GRID["w_air"]
        )
    )

    # Convert to list of tuples with descriptions
    configs = []
    for i, (vdw, elec, desolv, air) in enumerate(combinations):
        description = f"comb_{i:03d}_v{vdw}_e{elec}_d{desolv}_a{air}"
        configs.append((vdw, elec, desolv, air, description))

    return configs


def main():
    print("=" * 70)
    print("GDock Weight Optimization - COMBINATORIAL GRID SEARCH")
    print("=" * 70)

    # Generate all configurations
    configs = generate_configurations()
    total_runs = len(configs)

    print("Grid configuration:")
    print(f"  VDW:    {GRID['w_vdw']}")
    print(f"  Elec:   {GRID['w_elec']}")
    print(f"  Desolv: {GRID['w_desolv']}")
    print(f"  AIR:    {GRID['w_air']}")
    print()
    print(f"Total combinations: {total_runs}")

    print()

    # Confirm
    response = input(f"This will run {total_runs} configurations. Continue? [y/N]: ")
    if response.lower() != "y":
        print("Aborted.")
        sys.exit(0)

    print()

    # Build gdock
    build()

    # Create output directory
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)

    # Results storage
    results = []

    print("=" * 70)
    print()

    # Track progress
    start_time = datetime.now()

    # Run each configuration
    for i, (w_vdw, w_elec, w_desolv, w_air, description) in enumerate(configs, 1):
        print(
            f"[{i}/{total_runs}] VDW={w_vdw} Elec={w_elec} Desolv={w_desolv} AIR={w_air}"
        )

        # Log file for this run
        log_file = log_dir / f"grid_{description}.log"

        # Run gdock
        success = run(w_vdw, w_elec, w_desolv, w_air, log_file)

        if not success:
            print("  → Run failed, skipping...")
            print()
            continue

        # Extract results
        best_dockq = extract_best_dockq(log_file)
        best_score, best_rmsd, best_fnat = extract_final_metrics(log_file)
        conv_gen = extract_convergence_gen(log_file)

        print(f"  → DockQ: {best_dockq:.3f}")

        # Store results
        results.append(
            {
                "w_vdw": w_vdw,
                "w_elec": w_elec,
                "w_desolv": w_desolv,
                "w_air": w_air,
                "best_dockq": best_dockq,
                "best_score": best_score if best_score else "N/A",
                "best_rmsd": best_rmsd if best_rmsd else "N/A",
                "best_fnat": best_fnat if best_fnat else "N/A",
                "convergence_gen": conv_gen if conv_gen else "N/A",
            }
        )

        # Save intermediate results every 10 runs
        if i % 10 == 0 or i == total_runs:
            df_temp = pd.DataFrame(results)
            df_temp.to_csv("weight_sweep_results.csv", index=False)
            print(f"  → Progress saved ({len(results)} runs completed)")

        print()

    # Final save
    df = pd.DataFrame(results)
    df.to_csv("weight_sweep_results.csv", index=False)

    total_time = datetime.now() - start_time
    hours = int(total_time.total_seconds() // 3600)
    minutes = int((total_time.total_seconds() % 3600) // 60)

    print("=" * 70)
    print("Combinatorial grid search completed!")
    print("=" * 70)
    print()
    print(f"Total time: {hours}h {minutes}m")
    print(f"Configurations tested: {len(results)}")
    print()

    # Show top 10 configurations
    print("Top 10 configurations by DockQ:")
    print("-" * 70)

    df_sorted = df.sort_values("best_dockq", ascending=False)
    print(f"{'Rank':<5} {'DockQ':<8} {'VDW':<6} {'Elec':<6} {'Desolv':<8} {'AIR':<8}")
    print("-" * 70)

    for idx, (_, row) in enumerate(df_sorted.head(10).iterrows(), 1):
        print(
            f"{idx:<5} {row['best_dockq']:<8.3f} {row['w_vdw']:<6.1f} "
            f"{row['w_elec']:<6.1f} {row['w_desolv']:<8.1f} {row['w_air']:<8.1f}"
        )

    print()
    print("Next step - analyze results:")
    print("  python3 analyze_sweep.py")
    print()

    # Save completion marker
    with open("sweep_completed.txt", "w") as f:
        f.write(f"Combinatorial grid search completed at: {datetime.now()}\n")
        f.write(f"Total runs: {len(results)}\n")
        f.write(f"Total time: {hours}h {minutes}m\n")
        f.write(f"Best DockQ: {df['best_dockq'].max():.3f}\n")
        f.write("\nGrid:\n")
        for key, values in GRID.items():
            f.write(f"  {key}: {values}\n")


if __name__ == "__main__":
    main()
