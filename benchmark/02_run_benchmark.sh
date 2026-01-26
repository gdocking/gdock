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

# Count total complexes for progress
total=$(ls -d "$DATA_DIR"/*/ 2>/dev/null | wc -l)
current=0

echo "Running benchmark on $total complexes"
if [[ -n "$NPROC" ]]; then
  echo "Using $NPROC processors per run"
else
  echo "Using default processors (total - 2)"
fi
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

  # Run gdock
  if "$GDOCK" run \
    --receptor "$receptor" \
    --ligand "$ligand" \
    --restraints "$restraints" \
    --reference "$reference" \
    --output-dir "$output_dir" \
    $NPROC_ARG \
    >"$output_dir.log" 2>&1; then

    # Extract DockQ from metrics.tsv
    if [[ -f "$output_dir/metrics.tsv" ]]; then
      dockq=$(awk -F'\t' 'NR==2 {print $2}' "$output_dir/metrics.tsv")
      echo "DockQ=$dockq"
    else
      echo "OK (no metrics)"
    fi
  else
    echo "FAILED (see $output_dir.log)"
  fi
done

echo ""
echo "Results saved to $RESULTS_DIR/"
