#!/bin/bash
# Example run script for gdock
# This example uses the provided test data (2OOB structure)

# Restraints are defined as specific pairs: receptor_residue:ligand_residue
# Each pair represents a distance restraint between those two residues
# These are actual interface residues from the native complex<D-d> (within 5.0Å)
# Chain A interface: 931, 933, 934, 936, 937, 938, 940, 941, 946, 950
# Chain B interface: 6, 8, 42, 44, 45, 46, 47, 48, 49, 66, 68, 69, 70

cargo run --release -- run \
  --receptor data/A.pdb \
  --ligand data/B.pdb \
  --restraints 933:6,936:8,940:42,941:44,946:45,950:46 \
  --reference data/2oob.pdb
