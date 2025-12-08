# Weight Calibration Pipeline

## Overview

This pipeline calibrates the energy function weights for gdock using the dockground bound docking decoys benchmark. The goal is to find the optimal combination of weights for VDW, electrostatics, and desolvation terms that best ranks near-native structures.

## tl:dr

Run the entire calibration pipeline:

```bash
# Step 0: Download dataset
./00_download_dataset.sh

# Step 1: Prepare structures
./01_prepare_dataset.sh

# Step 2: Extract raw energies
./02_extract_raw_energies.sh

# Step 3: Optimize weights
Rscript grid_search.R
```

## Dataset

**Source**: Dockground Model-Docking Decoy Set 1.0
**URL**: <https://dockground.compbio.ku.edu/downloads/model/decoy/model-docking-decoy-set1-1.0.tgz>

**Why this dataset?**:
This dataset is ideal for weight calibration because:

1. Each complex has exactly 1 near-native structure among 100 decoys
2. The scoring function should rank the near-native structure highest
3. We can optimize weights to maximize this ranking performance across all complexes

## Pipeline Steps

### Step 0: Download Dataset

Download the dockground decoy set (~2.4 GB):

```bash
./00_download_dataset.sh
```

**What it does**:

- Downloads `model-docking-decoy-set1-1.0.tgz` from dockground
- Extracts archive to `data/raw/`
- Verifies download and extraction

**Output**: `data/raw/` containing complex directories with decoy structures

### Step 1: Prepare Dataset

Extract and organize structures for scoring:

```bash
./01_prepare_dataset.sh
```

**Input**: `data/raw/{complex_directories}/`
**Output**: `data/prepared/{complex}/`

### Step 2: Extract Raw Energies

Score all complexes with weights = 1 to get raw energy values:

```bash
./02_extract_raw_energies.sh
```

- Runs gdock in `--score` mode with all weights set to 1
- No restraints (scoring pre-generated decoys, not guiding docking)
- Extracts VDW, electrostatics, and desolvation energies for each model

**Output**: `results/raw_energies.tsv`

Columns:

- `complex` - Complex identifier (e.g., 1a9n)
- `model` - Model number within ensemble
- `score` - Total score with equal weights (vdw + elec + desolv)
- `vdw` - Van der Waals energy
- `elec` - Electrostatic energy
- `desolv` - Desolvation energy
- `w_vdw`, `w_elec`, `w_desolv`, `w_air` - Weight values used (all 1)

### Step 3: Optimize Weights

Find the best weight combination for ranking models using grid search:

```bash
Rscript grid_search.R
```

1. Loads raw energies from `results/raw_energies.tsv`
2. Tests 174,050 weight combinations in parallel:
   - `w_vdw`: 0.1 to 5.0 (step 0.1) = 50 values
   - `w_elec`: 0.01 to 0.3 (step 0.005) = 59 values
   - `w_desolv`: 0.1 to 3.0 (step 0.05) = 59 values
3. For each combination:
   - Re-scores all models with new weights
   - Ranks models by score (lower is better)
   - Checks if model #1 (near-native) ranks in top-1, top-5, top-10, top-20, top-50
   - Calculates mean and median rank of model #1 across all complexes
4. Selects best weights by optimizing for:
   - Primary: top50 success rate (% of complexes where model #1 ranks in top 50)
   - Tiebreakers: top20, top10, top5, then mean_rank

**Output**:

- `results/grid_search_results.tsv` - Full results for all 174,050 combinations
- `results/optimized_weights.json` - Best weights with performance metrics

## Notes

- Dataset: Dockground Model-Docking Decoy Set 1.0 (100 decoys per complex)
- Near-native structure is always model #1 (acceptable or better by CAPRI criteria)
- AIR weights are not calibrated here since we don't use restraints when scoring pre-generated decoys
- Grid search uses parallel processing via R's `parallel` package (base R, no external dependencies)
- This is a completely reproducible pipeline from data download to final weight optimization
