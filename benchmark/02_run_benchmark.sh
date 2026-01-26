#!/bin/bash
# Run gdock on all benchmark complexes
#
# Usage: ./02_run_benchmark.sh [nproc]
#   nproc: number of processors to use (default: total - 2)

DATA_DIR="data"
RESULTS_DIR="results"
GDOCK="../target/release/gdock"

# Number of processors - override via command line or edit here for different machines
NPROC="${1:-}"

# Create results directory
mkdir -p "$RESULTS_DIR"

# Determine actual number of processors being used
if [[ -n "$NPROC" ]]; then
  NPROC_USED="$NPROC"
else
  total_cpus=$(nproc)
  NPROC_USED=$((total_cpus - 2))
  [[ $NPROC_USED -lt 1 ]] && NPROC_USED=1
fi

# Initialize timing summary file
echo -e "complex\tdockq\ttime_s\trec_atoms\tlig_atoms\trestraints\tnproc" > "$RESULTS_DIR/timing.tsv"

# Count total complexes for progress
total=$(ls -d "$DATA_DIR"/*/ 2>/dev/null | wc -l)
current=0

echo "Running benchmark on $total complexes"
echo "Using $NPROC_USED processors per run"
echo ""

for complex_dir in "$DATA_DIR"/*/; do
  pdb_id=$(basename "$complex_dir")
  ((current++))

  receptor="$complex_dir/receptor.pdb"
  ligand="$complex_dir/ligand.pdb"
  restraints="$complex_dir/restraints.txt"
  reference="$complex_dir/reference.pdb"
  output_dir="$RESULTS_DIR/$pdb_id"

  # Skip if missing files
  if [[ ! -f "$receptor" || ! -f "$ligand" || ! -f "$restraints" ]]; then
    echo "[$current/$total] $pdb_id: SKIPPED (missing input files)"
    continue
  fi

  # Create reference by combining receptor + ligand (if not exists)
  if [[ ! -f "$reference" ]]; then
    cat "$receptor" "$ligand" >"$reference"
  fi

  echo -n "[$current/$total] $pdb_id: "

  # Build nproc argument if set
  NPROC_ARG=""
  if [[ -n "$NPROC" ]]; then
    NPROC_ARG="--nproc $NPROC"
  fi

  # Count atoms and restraints for timing analysis
  n_rec_atoms=$(grep -c "^ATOM" "$receptor" 2>/dev/null || echo 0)
  n_lig_atoms=$(grep -c "^ATOM" "$ligand" 2>/dev/null || echo 0)
  # Count restraints (number of colon-separated pairs, e.g., "39:47,40:48" = 2 pairs)
  n_restraints=$(grep -o ':' "$restraints" 2>/dev/null | wc -l || echo 0)

  # Run gdock with timing
  start_time=$(date +%s.%N)
  if "$GDOCK" run \
    --receptor "$receptor" \
    --ligand "$ligand" \
    --restraints "$restraints" \
    --reference "$reference" \
    --output-dir "$output_dir" \
    $NPROC_ARG \
    >"$output_dir.log" 2>&1; then
    end_time=$(date +%s.%N)
    elapsed=$(echo "$end_time - $start_time" | bc)

    # Extract DockQ from metrics.tsv (column 4: dockq)
    if [[ -f "$output_dir/metrics.tsv" ]]; then
      dockq=$(awk -F'\t' 'NR==2 {print $4}' "$output_dir/metrics.tsv")
      printf "DockQ=%.3f  time=%.1fs  (rec=%d lig=%d res=%d)\n" "$dockq" "$elapsed" "$n_rec_atoms" "$n_lig_atoms" "$n_restraints"
      # Append to timing summary
      echo -e "$pdb_id\t$dockq\t$elapsed\t$n_rec_atoms\t$n_lig_atoms\t$n_restraints\t$NPROC_USED" >> "$RESULTS_DIR/timing.tsv"
    else
      printf "OK  time=%.1fs\n" "$elapsed"
    fi
  else
    end_time=$(date +%s.%N)
    elapsed=$(echo "$end_time - $start_time" | bc)
    printf "FAILED  time=%.1fs (see $output_dir.log)\n" "$elapsed"
  fi
done

echo ""
echo "Results saved to $RESULTS_DIR/"
