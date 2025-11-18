#!/bin/bash
# Example run script for GDock

# This example uses the provided test data (2OOB structure)
# Receptor: Chain A (residues 929-972)
# Ligand: Chain B (residues 1-72)
# Reference: Native complex for validation

# Restraints are defined as specific pairs: receptor_residue:ligand_residue
# Each pair represents a distance restraint between those two residues

cargo build --release

./target/release/gdock run \
  --receptor data/A.pdb \
  --ligand data/B.pdb \
  --restraints 933:14,936:16,940:18,949:19,956:28,960:38,970:40 \
  --reference data/2oob.pdb \
  --w_vdw 1.0 \
  --w_elec 0.5 \
  --w_desolv 0.5 \
  --w_air 100.0

echo ""
echo "Output files:"
echo "  - best_by_score.pdb: Best solution by energy score"
echo "  - best_by_dockq.pdb: Best solution by DockQ metric"
