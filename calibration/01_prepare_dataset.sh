#!/bin/bash
# Prepare dockground 0A (bound) decoy structures for scoring
#
# This script extracts the receptor and ligand files from the dockground
# bound docking decoys (0A RMSD = native/bound structures).
#
# Dataset structure:
#   - data/raw/*_0A.tar.gz contains bound structures
#   - Each archive has: {complex}0A_0A.pdb (receptor) and {complex}0B_0A_poses.pdb (ligand ensemble)

set -euo pipefail

RAW_DIR="data/ModelModel"
OUTPUT_DIR="data/prepared"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if raw data exists
if [[ ! -d "$RAW_DIR" ]]; then
  echo "ERROR: Raw data directory not found: $RAW_DIR"
  echo "Run ./00_download_dataset.sh first"
  exit 1
fi

# Find all 0A (bound) archives
archives=$(ls "$RAW_DIR"/*_0A.tar.gz 2>/dev/null || true)

if [[ -z "$archives" ]]; then
  echo "ERROR: No *_0A.tar.gz archives found in $RAW_DIR"
  echo "Expected bound (0A) decoy archives"
  exit 1
fi

echo $archives

for archive in $archives; do

  # Extract complex name (e.g., 1a9n from 1a9n_0A.tar.gz)
  basename_file=$(basename "$archive" _0A.tar.gz)
  echo " Processing $basename_file..."

  # Create output directory for this complex
  output_complex_dir="$OUTPUT_DIR/$basename_file"
  mkdir -p "$output_complex_dir"
  tar -xzf "$archive" -C "$output_complex_dir" 2>/dev/null

done

echo
echo "Done! Prepared structures in: $OUTPUT_DIR"
echo "Total complexes: $(ls -1d "$OUTPUT_DIR"/*/ 2>/dev/null | wc -l)"
