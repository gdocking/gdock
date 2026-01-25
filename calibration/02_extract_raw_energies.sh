#!/bin/bash
# Extract raw energy values from all complexes
#
# This script scores all prepared complexes with weights = 1 to extract
# the raw energy component values (VDW, elec, desolv).
#
# Output: results/raw_energies.tsv

set -euo pipefail
shopt -s nullglob

PREPARED_DIR="data/prepared"
OUTPUT_DIR="results"
GDOCK="../target/release/gdock"

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
echo

# Initialize output file with header
echo -e "complex\tmodel\tscore\tvdw\telec\tdesolv\tclash_pct\tw_vdw\tw_elec\tw_desolv\tw_air" >"$output_file"

for complex_dir in "$PREPARED_DIR"/*/; do

  complex=$(basename "$complex_dir")
  receptor=$(ls "${complex_dir}${complex}_0A/"*_0A.pdb)
  ligand=$(ls "${complex_dir}${complex}_0A/"*_0A_poses.pdb)

  if [[ ! -f $receptor ]] || [[ ! -f $ligand ]]; then
    echo " WARNING COULD NOT EVALUTE ${complex}"
    continue
  fi

  echo " Scoring $complex..."

  # Score with all weights = 1 to get raw energy values
  "$GDOCK" score \
    --receptor "$receptor" \
    --ligand "$ligand" \
    --w_vdw 1 \
    --w_elec 1 \
    --w_desolv 1 \
    --w_air 1 |
    tail -n +2 |
    awk -v complex="$complex" '{print complex "\t" $0}' >>"$output_file"
done

echo
echo "Done! Raw energies extracted to: $output_file"
