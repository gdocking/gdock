use crate::structure;
use rand::distributions::{Alphanumeric, DistString};

pub fn generate_restraints(
    molecule1: &structure::Molecule,
    molecule2: &structure::Molecule,
) -> (Vec<i16>, Vec<i16>) {
    let mut restraints_a = Vec::new();
    let mut restraints_b = Vec::new();

    // Calculate the distances between all atoms of molecule a and b
    for atom1 in &molecule1.0 {
        for atom2 in &molecule2.0 {
            let dist = structure::distance(atom1, atom2);
            // println!("Distance between {} and {} is {}", atom1.resnum, atom2.resnum, dist);
            if dist < 5.0 {
                // let mut r = Restraint::new();
                // r.resseq = atom1.resseq;
                restraints_a.push(atom1.resseq);

                // let mut r = Restraint::new();
                // r.resseq = atom2.resseq;
                restraints_b.push(atom2.resseq);
            }
        }
    }

    // Remove duplicates
    restraints_a.sort();
    restraints_a.dedup();

    restraints_b.sort();
    restraints_b.dedup();

    println!("Restraints A: {:?}", restraints_a);
    println!("Restraints B: {:?}", restraints_b);

    (restraints_a, restraints_b)
}

// Generate a unique string
pub fn unique_str() -> String {
    Alphanumeric.sample_string(&mut rand::thread_rng(), 16)
}

// pub fn is_converging(arr: &Vec<f64>) -> bool {
//     // Check if the last 10 values are within 0.01 of the last value
//     let mut last = arr[arr.len() - 1];

//     last = (last * 100.0).round() / 100.0;
//     let mut count = 0;
//     for i in (arr.len() - 10)..arr.len() {
//         if (arr[i] - last).abs() < 0.1 {
//             count += 1;
//         }
//         last = arr[i];
//     }
//     count == 10
// }

pub fn position_ligand(
    receptor: &structure::Molecule,
    mut ligand: structure::Molecule,
) -> structure::Molecule {
    // Position the ligand in the center of the receptor
    let (r_x, r_y, r_z) = receptor.center_of_mass();
    let (l_x, l_y, l_z) = ligand.center_of_mass();

    // Move the ligand to the center of the receptor
    for atom in ligand.0.iter_mut() {
        atom.x -= l_x;
        atom.y -= l_y;
        atom.z -= l_z;
    }

    ligand.translate(r_x, r_y, r_z);

    ligand
}
