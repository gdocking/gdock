// use crate::structure::write_pdb
use crate::restraints;
use crate::structure;
use crate::structure::Atom;

// use crate::structure;

// // Calculate the root mean square deviation between two molecules
// pub fn calc_rmsd(molecule1: &structure::Molecule, molecule2: &structure::Molecule) -> f64 {
//     let n = molecule1.0.len();
//     assert_eq!(
//         n,
//         molecule2.0.len(),
//         "Molecules must have the same number of atoms"
//     );

//     // structure::write_pdb(&molecule1, &"model_rmsd.pdb".to_string());
//     // // utils::write_pdb(&self.receptor, &"receptor.pdb".to_string());
//     // structure::write_pdb(&molecule2, &"reference_rmsd.pdb".to_string());

//     // utils::write_pdb();

//     // Calculate the sum of squared distances between corresponding atoms
//     let sum_squared_distances = molecule1
//         .0
//         .iter()
//         .zip(molecule2.0.iter())
//         .map(|(atom1, atom2)| atom1.distance_to(atom2))
//         .sum::<f64>();

//     // Calculate the root mean square deviation

//     // println!("{}", rmsd);

//     // panic!();

//     (sum_squared_distances / n as f64).sqrt()
// }

// Calculate the ratio of clashes
//   i.e. the number of atom pairs that are closer than the sum of their van der Waals radii
pub fn calc_clashes(molecule1: &structure::Molecule, molecule2: &structure::Molecule) -> f64 {
    let mut clashes = 0;
    let mut total_atoms = 0;

    for atom1 in &molecule1.0 {
        for atom2 in &molecule2.0 {
            total_atoms += 1;

            let dist = structure::distance(atom1, atom2);
            let vdw_dist = atom1.vdw_radius + atom2.vdw_radius;
            if dist < vdw_dist {
                clashes += 1;
            }
        }
    }
    clashes as f64 / total_atoms as f64
}

/// Calculates the Lennard-Jones potential energy between two particles.
///
/// The Lennard-Jones potential energy is a pairwise interaction potential used to model
/// the intermolecular forces between particles. It is commonly used in molecular dynamics simulations.
fn lj_potential(atom1: &Atom, atom2: &Atom, distance: f64) -> f64 {
    let epsilon_ij = (atom1.epsilon * atom2.epsilon).sqrt();
    let rmin_ij = atom1.rmin2 + atom2.rmin2;
    let ratio_pow = (rmin_ij / distance).powi(6);
    epsilon_ij * (ratio_pow.powi(2) - 2.0 * ratio_pow)
}

/// Calculates the van der Waals energy between two molecules.
///
/// # Arguments
///
/// * `molecule1` - The first molecule.
/// * `molecule2` - The second molecule.
///
/// # Returns
///
/// The van der Waals energy between the two molecules.
pub fn vdw_energy(receptor: &structure::Molecule, ligand: &structure::Molecule) -> f64 {
    let mut energy = 0.0;

    for atom1 in &receptor.0 {
        for atom2 in &ligand.0 {
            // Check if atoms are not from the same molecule
            if atom1.serial != atom2.serial {
                let dist = structure::distance(atom1, atom2);
                if dist > 0.0 {
                    energy += lj_potential(atom1, atom2, dist);
                }
            }
        }
    }
    energy
}

// pub fn coulombic_energy(molecule1: &structure::Molecule, molecule2: &structure::Molecule) -> f64 {
//     let k = 8.9875517923e9; // Coulomb's constant in N m^2 C^-2
//     let mut energy = 0.0;

//     for atom1 in &molecule1.0 {
//         for atom2 in &molecule2.0 {
//             let dist = structure::distance(atom1, atom2);
//             if dist > 0.0 {
//                 let charge1 = atom1.charge.parse::<f64>().unwrap();
//                 let charge2 = atom2.charge.parse::<f64>().unwrap();
//                 energy += k * (charge1 * charge2) / dist;
//             }
//         }
//     }

//     energy
// }

// pub

// // Evaluate the ratio of restraints that are satisfied
// pub fn eval_restraints(molecule1: &structure::Molecule, molecule2: &structure::Molecule) -> f64 {
//     let cutoff = 10.0;
//     // Get the atoms which are part of a restraint
//     let mut restraints_a: Vec<&Atom> = Vec::new();
//     for atom in &molecule1.0 {
//         if molecule1.1.contains(&atom.resseq) && atom.name == "CA" {
//             restraints_a.push(atom);
//         }
//     }
//     let mut restraints_b: Vec<&Atom> = Vec::new();
//     for atom in &molecule2.0 {
//         if molecule2.1.contains(&atom.resseq) && atom.name == "CA" {
//             restraints_b.push(atom);
//         }
//     }

//     // let mut in_contact_mol_1: Vec<&i16> = Vec::new();

//     let mut in_contact_a = 0;

//     // println!("{} {}", restraints_a.len(), restraints_b.len());
//     // println!("{:?}", restraints_a);
//     // panic!();

//     for atom_a in &restraints_a {
//         for atom_b in &restraints_b {
//             let dist = structure::distance(atom_a, atom_b);
//             if dist <= cutoff {
//                 in_contact_a += 1;
//                 break;
//             }
//         }
//     }

//     let mut in_contact_b = 0;

//     for atom_b in &restraints_b {
//         for atom_a in &restraints_a {
//             let dist = structure::distance(atom_b, atom_a);
//             if dist <= cutoff {
//                 // println!("{} {}", atom_b.resseq, atom_a.resseq);
//                 // panic!("");
//                 in_contact_b += 1;
//                 break;
//             }
//         }
//     }

//     // Do the calculation
//     let v_a = in_contact_a as f64;
//     let v_b = in_contact_b as f64;

//     let t_a = molecule1.1.len() as f64;
//     let t_b = molecule2.1.len() as f64;

//     // println("{}")
//     // if score > 0.0 {
//     //     println!("{} {} {} {} {}", v_a, v_b, t_a, t_b, score);
//     //     panic!();
//     // }

//     (v_a + v_b) / (t_a + t_b)
// }

pub fn desolv_energy(_molecule1: &structure::Molecule, _molecule2: &structure::Molecule) -> f64 {
    todo!();
}

pub fn bsa_energy(_molecule1: &structure::Molecule, _molecule2: &structure::Molecule) -> f64 {
    todo!();
}

// pub fn air_energy(_molecule1: &structure::Molecule, _molecule2: &structure::Molecule) -> f64 {
//     todo!();
// }

// pub fn calculate_eair(molecule1: &structure::Molecule, molecule2: &structure::Molecule) -> f64 {
//     let active_atoms_a = molecule1
//         .0
//         .iter()
//         .filter(|atom| molecule1.1.contains(&atom.resseq))
//         .collect::<Vec<&Atom>>();

//     let active_atoms_b = molecule2
//         .0
//         .iter()
//         .filter(|atom| molecule2.1.contains(&atom.resseq))
//         .collect::<Vec<&Atom>>();

//     let mut sum_inverse_sixth_powers = 0.0;

//     for &atom_a in &active_atoms_a {
//         for &atom_b in &active_atoms_b {
//             let dist = structure::distance(atom_a, atom_b);
//             if dist > 0.0 && dist <= 2.0 {
//                 sum_inverse_sixth_powers += 1.0 / dist.powi(6);
//             }
//         }
//     }

//     if sum_inverse_sixth_powers > 0.0 {
//         sum_inverse_sixth_powers.powf(-1.0 / 6.0)
//     } else {
//         f64::INFINITY
//     }
// }

pub fn satisfaction_ratio(
    restraints: &[restraints::Restraint],
    receptor: &structure::Molecule,
    ligand: &structure::Molecule,
) -> f64 {
    let satisfied = restraints
        .iter()
        .filter(|&x| x.is_satisfied(receptor, ligand))
        .count();

    let total = restraints.len();
    if total > 0 {
        satisfied as f64 / total as f64
    } else {
        0.0
    }
}

pub fn coulombic_energy(receptor: &structure::Molecule, ligand: &structure::Molecule) -> f64 {
    let k = 8.9875517923e9; // Coulomb's constant in N m^2 C^-2
    let mut energy = 0.0;

    for atom1 in &receptor.0 {
        for atom2 in &ligand.0 {
            let dist = structure::distance(atom1, atom2);
            if dist > 0.0 {
                let charge1 = atom1.charge;
                let charge2 = atom2.charge;
                energy += k * (charge1 * charge2) / (dist * dist);
            }
        }
    }

    energy
}
