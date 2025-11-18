// use crate::structure;
// use std::collections::HashSet;

// pub struct Complex<'a> {
//     receptor: &'a structure::Molecule,
//     ligand: &'a structure::Molecule,
//     reference: &'a structure::Molecule,
// }

// #[derive(PartialEq, Debug, Hash, Eq)]
// struct Contact<'a> {
//     chain_i: char,
//     res_i: i16,
//     atom_i: &'a str,
//     chain_j: char,
//     res_j: i16,
//     atom_j: &'a str,
// }

// impl<'a> Complex<'_> {
//     pub fn new(
//         receptor: &'a structure::Molecule,
//         ligand: &'a structure::Molecule,
//         reference: &'a structure::Molecule,
//     ) -> Complex<'a> {
//         // todo!("Calculate the interface here!");
//         // Then re-use the interface to avoid calculating it again in fnat()
//         Complex {
//             receptor,
//             ligand,
//             reference,
//         }
//     }

//     pub fn rmsd(&self) -> f64 {
//         // Calculate the sum of squared distances between corresponding atoms
//         let n = self.reference.0.len();
//         let sum_squared_distances = &self
//             .ligand
//             .0
//             .iter()
//             .zip(self.reference.0.iter())
//             .map(|(atom1, atom2)| atom1.distance_to(atom2))
//             .sum::<f64>();
//         (sum_squared_distances / n as f64).sqrt()
//     }

//     pub fn fnat(&self) -> f64 {
//         // todo!("Optimize this function!");
//         // Calculate native contacts
//         let mut native_contacts: Vec<Contact> = Vec::new();
//         let mut docked_contacts: Vec<Contact> = Vec::new();

//         for atom in &self.receptor.0 {
//             for ref_atom in &self.reference.0 {
//                 let distance = atom.distance_to(ref_atom);
//                 if distance <= 5.0 {
//                     let c = Contact {
//                         chain_i: atom.chainid,
//                         res_i: atom.resseq,
//                         atom_i: &atom.name,
//                         chain_j: ref_atom.chainid,
//                         res_j: ref_atom.resseq,
//                         atom_j: &ref_atom.name,
//                     };
//                     native_contacts.push(c);
//                 }
//             }
//             for ref_atom in &self.ligand.0 {
//                 let distance = atom.distance_to(ref_atom);
//                 if distance <= 5.0 {
//                     let c = Contact {
//                         chain_i: atom.chainid,
//                         res_i: atom.resseq,
//                         atom_i: &atom.name,
//                         chain_j: ref_atom.chainid,
//                         res_j: ref_atom.resseq,
//                         atom_j: &ref_atom.name,
//                     };
//                     docked_contacts.push(c);
//                 }
//             }
//         }

//         // Find out the percentage of the docked contacts that are in the native contacts
//         // Assuming Contact implements Eq and Hash
//         let native_contacts_set: HashSet<_> = native_contacts.into_iter().collect();
//         let mut n_native_contacts = 0;

//         for docked_contact in &docked_contacts {
//             if native_contacts_set.contains(docked_contact) {
//                 n_native_contacts += 1;
//             }
//         }
//         // let native_contacts_set: HashSet<_> = native_contacts.into_iter().collect();
//         // let mut n_native_contacts = 0;
//         // for docked_contact in &docked_contacts {
//         //     for native_contact in &native_contacts {
//         //         if docked_contact == native_contact {
//         //             n_native_contacts += 1;
//         //         }
//         //     }
//         // }

//         let n_native_contacts = n_native_contacts as f64;
//         let n_docked_contacts = docked_contacts.len() as f64;

//         if n_native_contacts == 0.0 || n_docked_contacts == 0.0 {
//             return 0.0;
//         }
//         n_native_contacts / n_docked_contacts
//     }

//     // fn irmsd() {
//     //     todo!()
//     // }

//     // fn lrmsd() {
//     //     todo!()
//     // }
// }
