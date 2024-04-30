use crate::structure;

#[derive(Debug, Clone)]
pub struct Restraint(structure::Atom, structure::Atom);

impl Restraint {
    fn new(atom1: structure::Atom, atom2: structure::Atom) -> Self {
        Self(atom1, atom2)
    }

    pub fn is_satisfied(
        &self,
        receptor: &structure::Molecule,
        ligand: &structure::Molecule,
    ) -> bool {
        // Look into mol1 for the atoms with the same resseq as self.0
        let i = receptor
            .0
            .iter()
            .filter(|&x| x.resseq == self.0.resseq)
            .collect::<Vec<_>>();

        // i must not be 0
        assert!(i.len() > 1);

        let j = ligand
            .0
            .iter()
            .filter(|&x| x.resseq == self.1.resseq)
            .collect::<Vec<_>>();

        // i and j must b
        assert!(j.len() > 1);

        // Loop over the atoms in i and j and check if ANY of the distances are lower than 5A
        for atom1 in i.iter() {
            for atom2 in j.iter() {
                if structure::distance(atom1, atom2) < 5.0 {
                    return true;
                }
            }
        }
        false
    }
}

pub fn create_restraints(mol1: &structure::Molecule, mol2: &structure::Molecule) -> Vec<Restraint> {
    let mut restraints = Vec::new();

    // Loop over the atoms in mol1
    for atom1 in mol1.0.iter() {
        // Loop over the atoms in mol2
        for atom2 in mol2.0.iter() {
            // Check if the atoms are within 5A
            if structure::distance(atom1, atom2) < 5.0 {
                // Create a new restraint
                let restraint = Restraint::new(atom1.clone(), atom2.clone());

                // Add the restraint to the vector
                restraints.push(restraint);
            }
        }
    }

    restraints
}
