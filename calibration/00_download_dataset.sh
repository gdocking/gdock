#!/bin/bash
# Download dockground model-docking decoy set
#
# This dataset is used to calibrate energy function weights by scoring decoys
# and optimizing weights to rank near-native structures highest.

set -euo pipefail

DATA_DIR="data"
DOWNLOAD_URL="https://dockground.compbio.ku.edu/downloads/model/decoy/model-docking-decoy-set1-1.0.tgz"
ARCHIVE_NAME="model-docking-decoy-set1-1.0.tgz"

echo "================================================"
echo "Downloading Dockground Model-Docking Decoy Set"
echo "================================================"
echo
echo "Dataset: Model-Docking Decoy Set 1.0"
echo "URL: $DOWNLOAD_URL"
echo "Description: 99 incorrect + 1 near-native decoys per complex"
echo

# Create data directory
mkdir -p "$DATA_DIR"

# Download dataset
if [[ -f "$DATA_DIR/$ARCHIVE_NAME" ]]; then
  echo "Archive already exists: $DATA_DIR/$ARCHIVE_NAME"
  echo "Skipping download. Delete the file to re-download."
else
  echo "Downloading dataset..."
  echo "This may take a while..."
  wget -O "$DATA_DIR/$ARCHIVE_NAME" "$DOWNLOAD_URL"
  echo
  echo "Download complete!"
fi

# Verify download
if [[ -f "$DATA_DIR/$ARCHIVE_NAME" ]]; then
  file_size=$(du -h "$DATA_DIR/$ARCHIVE_NAME" | cut -f1)
  echo
  echo "Dataset downloaded successfully!"
  echo "Location: $DATA_DIR/$ARCHIVE_NAME"
  echo "Size: $file_size"
else
  echo
  echo "ERROR: Download failed!"
  exit 1
fi

# Extract dataset
echo
echo "Extracting dataset..."
tar -xzf "$DATA_DIR/$ARCHIVE_NAME" -C "$DATA_DIR"

# Find the extracted directory
extracted_dir=$(tar -tzf "$DATA_DIR/$ARCHIVE_NAME" | head -1 | cut -d/ -f1)

if [[ -d "$DATA_DIR/$extracted_dir" ]]; then
  # Rename to standard 'raw' directory for consistency
  if [[ "$extracted_dir" != "raw" ]]; then
    mv "$DATA_DIR/$extracted_dir" "$DATA_DIR/raw"
    extracted_dir="raw"
  fi

  # Count complexes
  num_complexes=$(ls -1d "$DATA_DIR/$extracted_dir"/*/ 2>/dev/null | wc -l)

  echo "Extraction complete!"
  echo "Location: $DATA_DIR/$extracted_dir/"
  echo "Complexes found: $num_complexes"
else
  echo "ERROR: Extraction failed!"
  exit 1
fi

echo
echo "Next step: ./01_prepare_dataset.sh"
