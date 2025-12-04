use crate::structure;

#[derive(Debug, Clone)]
pub struct Restraint(pub structure::Atom, pub structure::Atom);

impl Restraint {
    fn new(atom1: structure::Atom, atom2: structure::Atom) -> Self {
        Self(atom1, atom2)
    }

    pub fn is_satisfied(
        &self,
        receptor: &structure::Molecule,
        ligand: &structure::Molecule,
    ) -> bool {
        // The restraint is based on CA atoms, so check CA-CA distance
        // This is consistent with how restraints are created

        // Find CA atom in receptor with matching resseq
        let ca_receptor = receptor
            .0
            .iter()
            .find(|x| x.resseq == self.0.resseq && x.name.trim() == "CA");

        // Find CA atom in ligand with matching resseq
        let ca_ligand = ligand
            .0
            .iter()
            .find(|x| x.resseq == self.1.resseq && x.name.trim() == "CA");

        // Both CAs must exist
        if let (Some(ca1), Some(ca2)) = (ca_receptor, ca_ligand) {
            let dist = structure::distance(ca1, ca2);
            // Restraint is satisfied if within flat-bottom bounds [0.0, 7.0]
            // This matches the upper_bound in air_energy()
            return dist <= 7.0;
        }

        false
    }
}

/// Create restraints from user-specified residue pairs
pub fn create_restraints_from_pairs(
    mol1: &structure::Molecule,
    mol2: &structure::Molecule,
    pairs: &[(i32, i32)],
) -> Vec<Restraint> {
    let mut restraints = Vec::new();

    for (res1, res2) in pairs {
        // Find CA atom in mol1 with matching residue number
        let ca_atom_1 = mol1
            .0
            .iter()
            .find(|atom| atom.name.trim() == "CA" && atom.resseq as i32 == *res1);

        // Find CA atom in mol2 with matching residue number
        let ca_atom_2 = mol2
            .0
            .iter()
            .find(|atom| atom.name.trim() == "CA" && atom.resseq as i32 == *res2);

        // Both CAs must exist
        if let (Some(atom1), Some(atom2)) = (ca_atom_1, ca_atom_2) {
            let restraint = Restraint::new(atom1.clone(), atom2.clone());
            restraints.push(restraint);
        } else {
            eprintln!(
                "Warning: Restraint pair {}:{} not found in structures",
                res1, res2
            );
        }
    }

    restraints
}

/// Create restraints automatically from interface (for benchmarking/development)
pub fn create_restraints(mol1: &structure::Molecule, mol2: &structure::Molecule) -> Vec<Restraint> {
    let mut restraints = Vec::new();

    // Only create restraints for CA atoms (backbone) to reduce the number of restraints
    // This is more realistic for protein-protein docking
    let ca_atoms_1: Vec<&structure::Atom> = mol1
        .0
        .iter()
        .filter(|atom| atom.name.trim() == "CA")
        .collect();

    let ca_atoms_2: Vec<&structure::Atom> = mol2
        .0
        .iter()
        .filter(|atom| atom.name.trim() == "CA")
        .collect();

    // Loop over CA atoms in mol1
    for atom1 in ca_atoms_1.iter() {
        // Loop over CA atoms in mol2
        for atom2 in ca_atoms_2.iter() {
            // Create restraints for CA atoms within 7.0A (interface contacts)
            // This matches the upper_bound used in air_energy() flat-bottom potential
            let dist = structure::distance(atom1, atom2);
            if dist < 7.0 {
                // Create a new restraint
                let restraint = Restraint::new((*atom1).clone(), (*atom2).clone());

                // Add the restraint to the vector
                restraints.push(restraint);
            }
        }
    }

    restraints
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_atom(name: &str, resseq: i16, x: f64, y: f64, z: f64) -> structure::Atom {
        structure::Atom {
            serial: 1,
            name: name.to_string(),
            altloc: ' ',
            resname: "ALA".to_string(),
            chainid: 'A',
            resseq,
            icode: ' ',
            x,
            y,
            z,
            occupancy: 1.0,
            tempfactor: 0.0,
            element: "C".to_string(),
            charge: 0.0,
            vdw_radius: 1.7,
            epsilon: 0.0,
            rmin2: 0.0,
            eps_1_4: 0.0,
            rmin2_1_4: 0.0,
        }
    }

    #[test]
    fn test_restraint_new() {
        let atom1 = create_test_atom("CA", 1, 0.0, 0.0, 0.0);
        let atom2 = create_test_atom("CA", 10, 5.0, 0.0, 0.0);

        let restraint = Restraint::new(atom1.clone(), atom2.clone());

        assert_eq!(restraint.0.resseq, 1);
        assert_eq!(restraint.1.resseq, 10);
    }

    #[test]
    fn test_restraint_is_satisfied_within_bounds() {
        let atom1 = create_test_atom("CA", 1, 0.0, 0.0, 0.0);
        let atom2 = create_test_atom("CA", 10, 5.0, 0.0, 0.0);

        let restraint = Restraint::new(atom1.clone(), atom2.clone());

        // Create molecules with CA atoms within 7.0A
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom("CA", 1, 0.0, 0.0, 0.0));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 10, 5.0, 0.0, 0.0)); // Distance = 5.0

        assert!(restraint.is_satisfied(&receptor, &ligand));
    }

    #[test]
    fn test_restraint_is_satisfied_at_boundary() {
        let atom1 = create_test_atom("CA", 1, 0.0, 0.0, 0.0);
        let atom2 = create_test_atom("CA", 10, 7.0, 0.0, 0.0);

        let restraint = Restraint::new(atom1.clone(), atom2.clone());

        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom("CA", 1, 0.0, 0.0, 0.0));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 10, 7.0, 0.0, 0.0)); // Distance = 7.0 (at boundary)

        assert!(restraint.is_satisfied(&receptor, &ligand));
    }

    #[test]
    fn test_restraint_is_not_satisfied_beyond_bounds() {
        let atom1 = create_test_atom("CA", 1, 0.0, 0.0, 0.0);
        let atom2 = create_test_atom("CA", 10, 10.0, 0.0, 0.0);

        let restraint = Restraint::new(atom1.clone(), atom2.clone());

        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom("CA", 1, 0.0, 0.0, 0.0));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 10, 10.0, 0.0, 0.0)); // Distance = 10.0

        assert!(!restraint.is_satisfied(&receptor, &ligand));
    }

    #[test]
    fn test_restraint_is_not_satisfied_missing_receptor_ca() {
        let atom1 = create_test_atom("CA", 1, 0.0, 0.0, 0.0);
        let atom2 = create_test_atom("CA", 10, 5.0, 0.0, 0.0);

        let restraint = Restraint::new(atom1, atom2);

        // Receptor has wrong residue number
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom("CA", 99, 0.0, 0.0, 0.0));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 10, 5.0, 0.0, 0.0));

        assert!(!restraint.is_satisfied(&receptor, &ligand));
    }

    #[test]
    fn test_restraint_is_not_satisfied_missing_ligand_ca() {
        let atom1 = create_test_atom("CA", 1, 0.0, 0.0, 0.0);
        let atom2 = create_test_atom("CA", 10, 5.0, 0.0, 0.0);

        let restraint = Restraint::new(atom1, atom2);

        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom("CA", 1, 0.0, 0.0, 0.0));

        // Ligand has wrong residue number
        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 99, 5.0, 0.0, 0.0));

        assert!(!restraint.is_satisfied(&receptor, &ligand));
    }

    #[test]
    fn test_create_restraints_from_pairs_single_pair() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom("CA", 1, 0.0, 0.0, 0.0));
        receptor.0.push(create_test_atom("CB", 1, 0.5, 0.0, 0.0)); // Non-CA atom

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 10, 5.0, 0.0, 0.0));
        ligand.0.push(create_test_atom("CB", 10, 5.5, 0.0, 0.0)); // Non-CA atom

        let pairs = vec![(1, 10)];
        let restraints = create_restraints_from_pairs(&receptor, &ligand, &pairs);

        assert_eq!(restraints.len(), 1);
        assert_eq!(restraints[0].0.resseq, 1);
        assert_eq!(restraints[0].1.resseq, 10);
    }

    #[test]
    fn test_create_restraints_from_pairs_multiple_pairs() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom("CA", 1, 0.0, 0.0, 0.0));
        receptor.0.push(create_test_atom("CA", 2, 1.0, 0.0, 0.0));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 10, 5.0, 0.0, 0.0));
        ligand.0.push(create_test_atom("CA", 11, 6.0, 0.0, 0.0));

        let pairs = vec![(1, 10), (2, 11)];
        let restraints = create_restraints_from_pairs(&receptor, &ligand, &pairs);

        assert_eq!(restraints.len(), 2);
        assert_eq!(restraints[0].0.resseq, 1);
        assert_eq!(restraints[0].1.resseq, 10);
        assert_eq!(restraints[1].0.resseq, 2);
        assert_eq!(restraints[1].1.resseq, 11);
    }

    #[test]
    fn test_create_restraints_from_pairs_missing_residue() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom("CA", 1, 0.0, 0.0, 0.0));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 10, 5.0, 0.0, 0.0));

        // Request a pair where one residue doesn't exist
        let pairs = vec![(1, 10), (999, 10)];
        let restraints = create_restraints_from_pairs(&receptor, &ligand, &pairs);

        // Should only create 1 restraint (the valid one)
        assert_eq!(restraints.len(), 1);
        assert_eq!(restraints[0].0.resseq, 1);
    }

    #[test]
    fn test_create_restraints_from_pairs_no_ca_atoms() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom("CB", 1, 0.0, 0.0, 0.0)); // Not CA

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CB", 10, 5.0, 0.0, 0.0)); // Not CA

        let pairs = vec![(1, 10)];
        let restraints = create_restraints_from_pairs(&receptor, &ligand, &pairs);

        // Should create 0 restraints since no CA atoms exist
        assert_eq!(restraints.len(), 0);
    }

    #[test]
    fn test_create_restraints_within_distance() {
        let mut mol1 = structure::Molecule::new();
        mol1.0.push(create_test_atom("CA", 1, 0.0, 0.0, 0.0));
        mol1.0.push(create_test_atom("CA", 2, 100.0, 0.0, 0.0)); // Far away

        let mut mol2 = structure::Molecule::new();
        mol2.0.push(create_test_atom("CA", 10, 5.0, 0.0, 0.0)); // Within 7.0A of res 1
        mol2.0.push(create_test_atom("CA", 11, 200.0, 0.0, 0.0)); // Far away

        let restraints = create_restraints(&mol1, &mol2);

        // Should only create 1 restraint (1-10 pair within 7.0A)
        assert_eq!(restraints.len(), 1);
        assert_eq!(restraints[0].0.resseq, 1);
        assert_eq!(restraints[0].1.resseq, 10);
    }

    #[test]
    fn test_create_restraints_at_boundary() {
        let mut mol1 = structure::Molecule::new();
        mol1.0.push(create_test_atom("CA", 1, 0.0, 0.0, 0.0));

        let mut mol2 = structure::Molecule::new();
        mol2.0.push(create_test_atom("CA", 10, 6.99, 0.0, 0.0)); // Just under 7.0A

        let restraints = create_restraints(&mol1, &mol2);

        assert_eq!(restraints.len(), 1);
    }

    #[test]
    fn test_create_restraints_beyond_distance() {
        let mut mol1 = structure::Molecule::new();
        mol1.0.push(create_test_atom("CA", 1, 0.0, 0.0, 0.0));

        let mut mol2 = structure::Molecule::new();
        mol2.0.push(create_test_atom("CA", 10, 10.0, 0.0, 0.0)); // Beyond 7.0A

        let restraints = create_restraints(&mol1, &mol2);

        assert_eq!(restraints.len(), 0);
    }

    #[test]
    fn test_create_restraints_filters_non_ca() {
        let mut mol1 = structure::Molecule::new();
        mol1.0.push(create_test_atom("CA", 1, 0.0, 0.0, 0.0));
        mol1.0.push(create_test_atom("CB", 1, 0.0, 0.0, 0.0)); // Non-CA
        mol1.0.push(create_test_atom("N", 1, 0.0, 0.0, 0.0)); // Non-CA

        let mut mol2 = structure::Molecule::new();
        mol2.0.push(create_test_atom("CA", 10, 5.0, 0.0, 0.0));
        mol2.0.push(create_test_atom("CB", 10, 5.0, 0.0, 0.0)); // Non-CA

        let restraints = create_restraints(&mol1, &mol2);

        // Should only create 1 restraint (CA-CA pair)
        assert_eq!(restraints.len(), 1);
        assert_eq!(restraints[0].0.name.trim(), "CA");
        assert_eq!(restraints[0].1.name.trim(), "CA");
    }

    #[test]
    fn test_create_restraints_empty_molecules() {
        let mol1 = structure::Molecule::new();
        let mol2 = structure::Molecule::new();

        let restraints = create_restraints(&mol1, &mol2);

        assert_eq!(restraints.len(), 0);
    }
}
