/// The evaluator module contains the evaluator struct and its methods.
use crate::structure;
use std::collections::HashSet;

#[derive(Clone)]
pub struct Evaluator {
    // receptor: structure::Molecule,
    reference: structure::Molecule,
    interface: (HashSet<i16>, HashSet<i16>),
    reference_interface: structure::Molecule,
    receptor_interface: structure::Molecule,
    native_contacts: Vec<Contact>,
}

#[derive(PartialEq, Debug, Hash, Eq, Clone)]
pub struct Contact {
    chain_i: char,
    res_i: i16,
    atom_i: String,
    chain_j: char,
    res_j: i16,
    atom_j: String,
}

#[derive(Debug)]
pub struct Metrics {
    pub rmsd: f64,
    pub irmsd: f64,
    pub fnat: f64,
    pub dockq: f64,
}

impl Metrics {
    pub fn rank(&self) -> String {
        if self.dockq >= 0.8 {
            "High".to_string()
        } else if self.dockq >= 0.5 {
            "Medium".to_string()
        } else {
            "Low".to_string()
        }
    }
}

impl Evaluator {
    pub fn new(receptor: structure::Molecule, reference: structure::Molecule) -> Evaluator {
        let interface = calculate_interface(&receptor, &reference);

        let receptor_interface = structure::filter_by_resseq_vec(&receptor, &interface.0);
        let reference_interface = structure::filter_by_resseq_vec(&reference, &interface.1);

        let native_contacts = calculate_contacts(&receptor_interface, &reference_interface);

        Evaluator {
            // receptor,
            reference,
            interface,
            reference_interface,
            receptor_interface,
            native_contacts,
        }
    }

    pub fn calc_metrics(&self, model: &structure::Molecule) -> Metrics {
        let fnat = self.calc_fnat(model);
        let irmsd = self.calc_irmsd(model);
        let rmsd = self.calc_rmsd(model);

        // Calculate DockQ using the same values
        let fnat_score = fnat;
        let irmsd_score = 1.0 / (1.0 + (irmsd / 1.5).powi(2));
        let lrmsd_score = 1.0 / (1.0 + (rmsd / 8.5).powi(2));
        let dockq = (fnat_score + irmsd_score + lrmsd_score) / 3.0;

        Metrics {
            rmsd,
            irmsd,
            fnat,
            dockq,
        }
    }

    fn calc_rmsd(&self, ligand: &structure::Molecule) -> f64 {
        // The receptor does not move, so we don't need to calculate over it,
        // then the RMSD here is the L-RMSD

        // Calculate the sum of squared distances between corresponding atoms
        let n = self.reference.0.len();
        let sum_squared_distances = &ligand
            .0
            .iter()
            .zip(self.reference.0.iter())
            .map(|(atom1, atom2)| atom1.distance_to(atom2))
            .sum::<f64>();
        (sum_squared_distances / n as f64).sqrt()
    }

    fn calc_irmsd(&self, ligand: &structure::Molecule) -> f64 {
        // The receptor does not move, so we don't need to calculate over it,
        // here we only need to calculate the rmsd between the reference and the ligand,
        // but only for the atoms that are in the interface.

        let target_ligand = structure::filter_by_resseq_vec(ligand, &self.interface.1);

        // Calculate the sum of squared distances between corresponding atoms
        let n = self.reference_interface.0.len();
        let sum_squared_distances = &target_ligand
            .0
            .iter()
            .zip(self.reference_interface.0.iter())
            .map(|(atom1, atom2)| atom1.distance_to(atom2))
            .sum::<f64>();
        (sum_squared_distances / n as f64).sqrt()
    }

    fn calc_fnat(&self, ligand: &structure::Molecule) -> f64 {
        // Filter ligand to only interface residues to match how native_contacts was calculated
        let ligand_interface = structure::filter_by_resseq_vec(ligand, &self.interface.1);

        // Get the contact between this ligand interface and the receptor interface
        let docked_contacts = calculate_contacts(&self.receptor_interface, &ligand_interface);

        // Check how many of the docked contacts are in the native contacts, this will be the common contacts
        let common_contacts = docked_contacts
            .iter()
            .filter(|c| self.native_contacts.contains(c))
            .count();

        // Calculate the FNAT
        common_contacts as f64 / self.native_contacts.len() as f64
    }
}

fn calculate_contacts(
    receptor: &structure::Molecule,
    ligand: &structure::Molecule,
) -> Vec<Contact> {
    let mut native_contacts: Vec<Contact> = Vec::new();
    for atom in &receptor.0 {
        for ref_atom in &ligand.0 {
            let distance = atom.distance_to(ref_atom);
            if distance <= 5.0 {
                let c = Contact {
                    chain_i: atom.chainid,
                    res_i: atom.resseq,
                    atom_i: atom.name.to_string(),
                    chain_j: ref_atom.chainid,
                    res_j: ref_atom.resseq,
                    atom_j: ref_atom.name.to_string(),
                };
                native_contacts.push(c);
            }
        }
    }
    native_contacts
}
// Check which residues are in contact between the receptor and the reference
fn calculate_interface(
    receptor: &structure::Molecule,
    reference: &structure::Molecule,
) -> (HashSet<i16>, HashSet<i16>) {
    let mut receptor_interface: HashSet<i16> = HashSet::new();
    let mut reference_interface: HashSet<i16> = HashSet::new();
    for atom in &receptor.0 {
        for ref_atom in &reference.0 {
            let distance = atom.distance_to(ref_atom);
            if distance <= 5.0 {
                receptor_interface.insert(atom.resseq);
                reference_interface.insert(ref_atom.resseq);
            }
        }
    }
    (receptor_interface, reference_interface)
}
