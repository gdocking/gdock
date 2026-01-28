use crate::structure::{read_pdb, Molecule};
use std::collections::HashSet;

/// Find interface residue pairs between receptor and ligand within a distance cutoff.
/// Returns pairs of (receptor_resseq, ligand_resseq) that are in contact.
pub fn find_interface_pairs(
    receptor: &Molecule,
    ligand: &Molecule,
    cutoff: f64,
) -> Vec<(i16, i16)> {
    let mut interface_pairs: Vec<(i16, i16)> = Vec::new();

    // Get unique residues from receptor and ligand
    let rec_residues: HashSet<i16> = receptor.0.iter().map(|a| a.resseq).collect();
    let lig_residues: HashSet<i16> = ligand.0.iter().map(|a| a.resseq).collect();

    for rec_res in &rec_residues {
        for lig_res in &lig_residues {
            // Get heavy atoms for this residue pair
            let rec_atoms: Vec<_> = receptor
                .0
                .iter()
                .filter(|a| a.resseq == *rec_res && !a.name.starts_with('H'))
                .collect();
            let lig_atoms: Vec<_> = ligand
                .0
                .iter()
                .filter(|a| a.resseq == *lig_res && !a.name.starts_with('H'))
                .collect();

            // Check if any pair of atoms is within cutoff
            'outer: for rec_atom in &rec_atoms {
                for lig_atom in &lig_atoms {
                    let dx = rec_atom.x - lig_atom.x;
                    let dy = rec_atom.y - lig_atom.y;
                    let dz = rec_atom.z - lig_atom.z;
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                    if dist <= cutoff {
                        interface_pairs.push((*rec_res, *lig_res));
                        break 'outer;
                    }
                }
            }
        }
    }

    // Sort for consistent output
    interface_pairs.sort();
    interface_pairs
}

/// Generate restraints from interface contacts and print to stdout
pub fn generate_restraints(receptor_file: String, ligand_file: String, cutoff: f64) {
    let receptor_model = read_pdb(&receptor_file);
    let ligand_model = read_pdb(&ligand_file);

    let receptor = &receptor_model.0[0];
    let ligand = &ligand_model.0[0];

    let interface_pairs = find_interface_pairs(receptor, ligand, cutoff);

    // Output in gdock restraints format
    let restraints_str: Vec<String> = interface_pairs
        .iter()
        .map(|(r, l)| format!("{}:{}", r, l))
        .collect();

    println!("{}", restraints_str.join(","));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_interface_pairs_with_test_data() {
        let receptor_model = read_pdb(&"data/2oob_A.pdb".to_string());
        let ligand_model = read_pdb(&"data/2oob_B.pdb".to_string());

        let receptor = &receptor_model.0[0];
        let ligand = &ligand_model.0[0];

        let pairs = find_interface_pairs(receptor, ligand, 5.0);

        // Should find some interface pairs
        assert!(!pairs.is_empty(), "Should find interface pairs");

        // All pairs should be sorted
        let mut sorted = pairs.clone();
        sorted.sort();
        assert_eq!(pairs, sorted, "Pairs should be sorted");
    }

    #[test]
    fn test_find_interface_pairs_cutoff_effect() {
        let receptor_model = read_pdb(&"data/2oob_A.pdb".to_string());
        let ligand_model = read_pdb(&"data/2oob_B.pdb".to_string());

        let receptor = &receptor_model.0[0];
        let ligand = &ligand_model.0[0];

        let pairs_5a = find_interface_pairs(receptor, ligand, 5.0);
        let pairs_10a = find_interface_pairs(receptor, ligand, 10.0);

        // Larger cutoff should find more or equal pairs
        assert!(
            pairs_10a.len() >= pairs_5a.len(),
            "10Å cutoff should find >= pairs than 5Å"
        );
    }
}
