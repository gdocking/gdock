#!/bin/bash
# Extract raw energy values from all complexes (PARALLEL VERSION)
#
# This script scores all prepared complexes with weights = 1 to extract
# the raw energy component values (VDW, elec, desolv, clash_pct).
#
# Output: results/raw_energies.tsv

set -euo pipefail
shopt -s nullglob

PREPARED_DIR="data/prepared"
OUTPUT_DIR="results"
GDOCK="../target/release/gdock"
TEMP_DIR=$(mktemp -d)
N_JOBS=$(( $(nproc) - 2 ))
# Ensure at least 1 job
[[ $N_JOBS -lt 1 ]] && N_JOBS=1

# Cleanup temp dir on exit
trap "rm -rf $TEMP_DIR" EXIT

# Build gdock if not already built
if [[ ! -f "$GDOCK" ]]; then
  echo "Building gdock..."
  cd ..
  cargo build --release
  cd calibration
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Output file
output_file="$OUTPUT_DIR/raw_energies.tsv"

echo "Output: $output_file"
echo "Using $N_JOBS parallel jobs"
echo

# Function to score a single complex
score_complex() {
  local complex_dir="$1"
  local gdock="$2"
  local temp_dir="$3"

  local complex=$(basename "$complex_dir")
  local receptor=$(ls "${complex_dir}${complex}_0A/"*_0A.pdb 2>/dev/null || true)
  local ligand=$(ls "${complex_dir}${complex}_0A/"*_0A_poses.pdb 2>/dev/null || true)

  if [[ ! -f "$receptor" ]] || [[ ! -f "$ligand" ]]; then
    echo "WARNING: Could not evaluate ${complex}" >&2
    return 0
  fi

  # Score with all weights = 1 to get raw energy values
  "$gdock" score \
    --receptor "$receptor" \
    --ligand "$ligand" \
    --w_vdw 1 \
    --w_elec 1 \
    --w_desolv 1 \
    --w_air 1 2>/dev/null |
    tail -n +2 |
    awk -v complex="$complex" '{print complex "\t" $0}' > "${temp_dir}/${complex}.tsv"

  echo "  ✓ $complex"
}

export -f score_complex

# Get list of complex directories
complex_dirs=("$PREPARED_DIR"/*/)

echo "Scoring ${#complex_dirs[@]} complexes..."
echo

# Run in parallel using GNU parallel if available, otherwise xargs
if command -v parallel &> /dev/null; then
  printf '%s\n' "${complex_dirs[@]}" | \
    parallel -j "$N_JOBS" score_complex {} "$GDOCK" "$TEMP_DIR"
else
  # Fallback to xargs (slightly less elegant but works)
  printf '%s\n' "${complex_dirs[@]}" | \
    xargs -P "$N_JOBS" -I {} bash -c 'score_complex "$@"' _ {} "$GDOCK" "$TEMP_DIR"
fi

echo
echo "Merging results..."

# Initialize output file with header
echo -e "complex\tmodel\tscore\tvdw\telec\tdesolv\tw_vdw\tw_elec\tw_desolv\tw_air" > "$output_file"

# Merge all temp files (sorted by complex name for reproducibility)
for f in $(ls "$TEMP_DIR"/*.tsv 2>/dev/null | sort); do
  cat "$f" >> "$output_file"
done

n_complexes=$(ls "$TEMP_DIR"/*.tsv 2>/dev/null | wc -l)
echo
echo "Done! Scored $n_complexes complexes"
echo "Results: $output_file"
