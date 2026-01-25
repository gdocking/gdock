#!/bin/bash
# Run gdock on all benchmark complexes

DATA_DIR="data"
RESULTS_DIR="results"
GDOCK="../target/release/gdock"

# Create results directory
mkdir -p "$RESULTS_DIR"

# Count total complexes for progress
total=$(ls -d "$DATA_DIR"/*/ 2>/dev/null | wc -l)
current=0

echo "Running benchmark on $total complexes"
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
        cat "$receptor" "$ligand" > "$reference"
    fi

    echo -n "[$current/$total] $pdb_id: "

    # Run gdock
    if "$GDOCK" run \
        --receptor "$receptor" \
        --ligand "$ligand" \
        --restraints "$restraints" \
        --reference "$reference" \
        --output-dir "$output_dir" \
        > "$output_dir.log" 2>&1; then

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
