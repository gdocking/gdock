/// The evaluator module contains the evaluator struct and its methods.
use crate::structure;
use std::collections::HashSet;

#[derive(Debug, Clone)]
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

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_atom(x: f64, y: f64, z: f64, resseq: i16, chainid: char) -> structure::Atom {
        structure::Atom {
            serial: 1,
            name: "CA".to_string(),
            altloc: ' ',
            resname: "ALA".to_string(),
            chainid,
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
            epsilon: -0.1,
            rmin2: 2.0,
            eps_1_4: -0.1,
            rmin2_1_4: 1.9,
        }
    }

    #[test]
    fn test_contact_equality() {
        let contact1 = Contact {
            chain_i: 'A',
            res_i: 1,
            atom_i: "CA".to_string(),
            chain_j: 'B',
            res_j: 10,
            atom_j: "CA".to_string(),
        };

        let contact2 = Contact {
            chain_i: 'A',
            res_i: 1,
            atom_i: "CA".to_string(),
            chain_j: 'B',
            res_j: 10,
            atom_j: "CA".to_string(),
        };

        assert_eq!(contact1, contact2);
    }

    #[test]
    fn test_metrics_rank_high() {
        let metrics = Metrics {
            rmsd: 1.0,
            irmsd: 0.5,
            fnat: 0.9,
            dockq: 0.85,
        };

        assert_eq!(metrics.rank(), "High");
    }

    #[test]
    fn test_metrics_rank_medium() {
        let metrics = Metrics {
            rmsd: 3.0,
            irmsd: 2.0,
            fnat: 0.6,
            dockq: 0.6,
        };

        assert_eq!(metrics.rank(), "Medium");
    }

    #[test]
    fn test_metrics_rank_low() {
        let metrics = Metrics {
            rmsd: 10.0,
            irmsd: 8.0,
            fnat: 0.2,
            dockq: 0.3,
        };

        assert_eq!(metrics.rank(), "Low");
    }

    #[test]
    fn test_calculate_contacts_within_distance() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 1, 'A'));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom(3.0, 0.0, 0.0, 10, 'B')); // Within 5.0Å

        let contacts = calculate_contacts(&receptor, &ligand);

        assert_eq!(contacts.len(), 1);
        assert_eq!(contacts[0].res_i, 1);
        assert_eq!(contacts[0].res_j, 10);
    }

    #[test]
    fn test_calculate_contacts_beyond_distance() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 1, 'A'));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom(10.0, 0.0, 0.0, 10, 'B')); // Beyond 5.0Å

        let contacts = calculate_contacts(&receptor, &ligand);

        assert_eq!(contacts.len(), 0);
    }

    #[test]
    fn test_calculate_interface() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 1, 'A'));
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 2, 'A'));
        receptor.0.push(create_test_atom(20.0, 0.0, 0.0, 3, 'A')); // Far away

        let mut reference = structure::Molecule::new();
        reference.0.push(create_test_atom(3.0, 0.0, 0.0, 10, 'B')); // Close to res 1 and 2

        let (receptor_interface, reference_interface) = calculate_interface(&receptor, &reference);

        assert_eq!(receptor_interface.len(), 2); // Res 1 and 2
        assert!(receptor_interface.contains(&1));
        assert!(receptor_interface.contains(&2));
        assert!(!receptor_interface.contains(&3));

        assert_eq!(reference_interface.len(), 1); // Res 10
        assert!(reference_interface.contains(&10));
    }

    #[test]
    fn test_evaluator_creation() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 1, 'A'));

        let mut reference = structure::Molecule::new();
        reference.0.push(create_test_atom(3.0, 0.0, 0.0, 10, 'B'));

        let evaluator = Evaluator::new(receptor, reference);

        assert_eq!(evaluator.interface.0.len(), 1); // Receptor has 1 interface residue
        assert_eq!(evaluator.interface.1.len(), 1); // Reference has 1 interface residue
    }

    #[test]
    fn test_calc_rmsd_identical() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 1, 'A'));

        let mut reference = structure::Molecule::new();
        reference.0.push(create_test_atom(3.0, 0.0, 0.0, 10, 'B'));

        let evaluator = Evaluator::new(receptor, reference.clone());

        // Same as reference
        let rmsd = evaluator.calc_rmsd(&reference);
        assert_eq!(rmsd, 0.0, "RMSD should be 0 for identical structures");
    }

    #[test]
    fn test_calc_rmsd_different() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 1, 'A'));

        let mut reference = structure::Molecule::new();
        reference.0.push(create_test_atom(3.0, 0.0, 0.0, 10, 'B'));

        let evaluator = Evaluator::new(receptor, reference);

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom(6.0, 0.0, 0.0, 10, 'B')); // 3Å away from reference

        let rmsd = evaluator.calc_rmsd(&ligand);
        // RMSD calculation: sqrt(sum_of_distances / n)
        // For single atom 3Å apart: sqrt(3.0 / 1) = sqrt(3.0) ≈ 1.732
        assert!(
            (rmsd - 3.0_f64.sqrt()).abs() < 1e-10,
            "RMSD should be sqrt(3.0)"
        );
    }

    #[test]
    fn test_calc_fnat_perfect() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 1, 'A'));

        let mut reference = structure::Molecule::new();
        reference.0.push(create_test_atom(3.0, 0.0, 0.0, 10, 'B'));

        let evaluator = Evaluator::new(receptor, reference.clone());

        // Identical to reference - should have 100% native contacts
        let fnat = evaluator.calc_fnat(&reference);
        assert_eq!(fnat, 1.0, "FNAT should be 1.0 for identical structures");
    }

    #[test]
    fn test_calc_fnat_no_contacts() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 1, 'A'));

        let mut reference = structure::Molecule::new();
        reference.0.push(create_test_atom(3.0, 0.0, 0.0, 10, 'B'));

        let evaluator = Evaluator::new(receptor, reference);

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom(20.0, 0.0, 0.0, 10, 'B')); // Far away - no contacts

        let fnat = evaluator.calc_fnat(&ligand);
        assert_eq!(fnat, 0.0, "FNAT should be 0.0 for no contacts");
    }

    #[test]
    fn test_calc_metrics_complete() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 1, 'A'));

        let mut reference = structure::Molecule::new();
        reference.0.push(create_test_atom(3.0, 0.0, 0.0, 10, 'B'));

        let evaluator = Evaluator::new(receptor, reference.clone());

        let metrics = evaluator.calc_metrics(&reference);

        assert_eq!(metrics.rmsd, 0.0);
        assert_eq!(metrics.irmsd, 0.0);
        assert_eq!(metrics.fnat, 1.0);
        assert!(
            metrics.dockq > 0.9,
            "DockQ should be high for native structure"
        );
    }

    #[test]
    fn test_dockq_calculation() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 1, 'A'));

        let mut reference = structure::Molecule::new();
        reference.0.push(create_test_atom(3.0, 0.0, 0.0, 10, 'B'));

        let evaluator = Evaluator::new(receptor, reference.clone());

        let metrics = evaluator.calc_metrics(&reference);

        // DockQ = (fnat + 1/(1+(irmsd/1.5)^2) + 1/(1+(lrmsd/8.5)^2)) / 3
        // For perfect match: (1 + 1 + 1) / 3 = 1.0
        assert!(
            (metrics.dockq - 1.0).abs() < 1e-10,
            "DockQ should be 1.0 for perfect match"
        );
    }
}
