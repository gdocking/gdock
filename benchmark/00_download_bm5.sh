#!/bin/bash
# Download and organize BM5.5 dataset

DOWNLOAD_URL="https://zlab.wenglab.org/benchmark/benchmark5.5.tgz"
ARCHIVE="benchmark5.5.tgz"
DATA_DIR="data"

# Download
if [[ ! -f "$ARCHIVE" ]]; then
  echo "Downloading BM5.5..."
  wget "$DOWNLOAD_URL"
else
  echo "Archive already exists, skipping download"
fi

# Extract
echo "Extracting..."
tar -xzf "$ARCHIVE"

# Create data directory
mkdir -p "$DATA_DIR"

# Reorganize: create per-complex folders with bound structures only
echo "Organizing structures..."
for receptor in benchmark5.5/structures/*_r_b.pdb; do
  # Extract PDB ID (e.g., 1AHW from 1AHW_r_b.pdb)
  pdb_id=$(basename "$receptor" _r_b.pdb)
  ligand="benchmark5.5/structures/${pdb_id}_l_b.pdb"

  if [[ -f "$ligand" ]]; then
    mkdir -p "$DATA_DIR/$pdb_id"
    cp "$receptor" "$DATA_DIR/$pdb_id/receptor.pdb"
    cp "$ligand" "$DATA_DIR/$pdb_id/ligand.pdb"
    echo "  $pdb_id"
  else
    echo "  WARNING: Missing ligand for $pdb_id"
  fi
done

# Create list of complexes
ls -1 "$DATA_DIR" >complexes.txt

# Summary
n_complexes=$(wc -l <complexes.txt)
echo ""
echo "Done! $n_complexes complexes organized in $DATA_DIR/"
echo "Complex list saved to complexes.txt"
