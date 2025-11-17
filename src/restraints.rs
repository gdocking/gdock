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
            // Use same threshold as restraint creation for consistency
            return dist < 8.0;
        }

        false
    }
}

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
            // Create restraints for CA atoms within 8A (native contacts)
            // This matches the typical interface definition
            let dist = structure::distance(atom1, atom2);
            if dist < 8.0 {
                // Create a new restraint
                let restraint = Restraint::new((*atom1).clone(), (*atom2).clone());

                // Add the restraint to the vector
                restraints.push(restraint);
            }
        }
    }

    restraints
}
