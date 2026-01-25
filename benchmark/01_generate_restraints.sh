#!/bin/bash
# Generate interface restraints for all complexes

DATA_DIR="data"
GDOCK="../target/release/gdock"
CUTOFF="${1:-5.0}"  # Default 5.0A, can be overridden via argument

echo "Generating restraints with cutoff=${CUTOFF}A"
echo ""

for complex_dir in "$DATA_DIR"/*/; do
    pdb_id=$(basename "$complex_dir")
    receptor="$complex_dir/receptor.pdb"
    ligand="$complex_dir/ligand.pdb"
    output="$complex_dir/restraints.txt"

    if [[ -f "$receptor" && -f "$ligand" ]]; then
        restraints=$("$GDOCK" restraints --receptor "$receptor" --ligand "$ligand" --cutoff "$CUTOFF")
        echo "$restraints" > "$output"
        n_pairs=$(echo "$restraints" | tr ',' '\n' | wc -l)
        echo "  $pdb_id: $n_pairs pairs"
    else
        echo "  $pdb_id: SKIPPED (missing files)"
    fi
done

echo ""
echo "Done! Restraints saved to {complex}/restraints.txt"
