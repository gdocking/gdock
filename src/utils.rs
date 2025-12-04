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

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_atom(x: f64, y: f64, z: f64) -> structure::Atom {
        structure::Atom {
            serial: 1,
            name: "CA".to_string(),
            altloc: ' ',
            resname: "ALA".to_string(),
            chainid: 'A',
            resseq: 1,
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
    fn test_unique_str_length() {
        let s = unique_str();
        assert_eq!(s.len(), 16);
    }

    #[test]
    fn test_unique_str_unique() {
        let s1 = unique_str();
        let s2 = unique_str();
        // With very high probability, two random 16-char strings will be different
        assert_ne!(s1, s2);
    }

    #[test]
    fn test_position_ligand_centers_overlap() {
        // Create receptor centered at (5, 10, 15)
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(3.0, 8.0, 13.0));
        receptor.0.push(create_test_atom(7.0, 12.0, 17.0));

        // Create ligand centered at (0, 0, 0)
        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom(-1.0, -1.0, -1.0));
        ligand.0.push(create_test_atom(1.0, 1.0, 1.0));

        let positioned = position_ligand(&receptor, ligand);

        // After positioning, ligand center should match receptor center
        let (cx, cy, cz) = positioned.center_of_mass();
        let (rx, ry, rz) = receptor.center_of_mass();

        assert!((cx - rx).abs() < 1e-10);
        assert!((cy - ry).abs() < 1e-10);
        assert!((cz - rz).abs() < 1e-10);
    }

    #[test]
    fn test_generate_restraints_within_distance() {
        let mut mol_a = structure::Molecule::new();
        let mut atom1 = create_test_atom(0.0, 0.0, 0.0);
        atom1.resseq = 1;
        mol_a.0.push(atom1);

        let mut mol_b = structure::Molecule::new();
        let mut atom2 = create_test_atom(3.0, 0.0, 0.0); // Distance = 3.0, within 5.0
        atom2.resseq = 10;
        mol_b.0.push(atom2);

        let (rest_a, rest_b) = generate_restraints(&mol_a, &mol_b);

        assert_eq!(rest_a.len(), 1);
        assert_eq!(rest_b.len(), 1);
        assert_eq!(rest_a[0], 1);
        assert_eq!(rest_b[0], 10);
    }

    #[test]
    fn test_generate_restraints_beyond_distance() {
        let mut mol_a = structure::Molecule::new();
        let mut atom1 = create_test_atom(0.0, 0.0, 0.0);
        atom1.resseq = 1;
        mol_a.0.push(atom1);

        let mut mol_b = structure::Molecule::new();
        let mut atom2 = create_test_atom(10.0, 0.0, 0.0); // Distance = 10.0, beyond 5.0
        atom2.resseq = 10;
        mol_b.0.push(atom2);

        let (rest_a, rest_b) = generate_restraints(&mol_a, &mol_b);

        assert_eq!(rest_a.len(), 0);
        assert_eq!(rest_b.len(), 0);
    }

    #[test]
    fn test_generate_restraints_deduplication() {
        let mut mol_a = structure::Molecule::new();
        let mut atom1 = create_test_atom(0.0, 0.0, 0.0);
        atom1.resseq = 1;
        let mut atom2 = create_test_atom(0.1, 0.0, 0.0);
        atom2.resseq = 1; // Same residue
        mol_a.0.push(atom1);
        mol_a.0.push(atom2);

        let mut mol_b = structure::Molecule::new();
        let mut atom3 = create_test_atom(1.0, 0.0, 0.0);
        atom3.resseq = 10;
        mol_b.0.push(atom3);

        let (rest_a, rest_b) = generate_restraints(&mol_a, &mol_b);

        // Should have only one entry even though two atoms from resseq 1 are close
        assert_eq!(rest_a.len(), 1);
        assert_eq!(rest_b.len(), 1);
    }
}
