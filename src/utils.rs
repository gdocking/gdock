use crate::structure;
use rand::distributions::{Alphanumeric, DistString};
use std::path::PathBuf;

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

/// Creates a unique temporary directory path.
///
/// This function generates a unique temporary directory path using the system's
/// temporary directory, the process ID, and a high-resolution timestamp. The directory
/// is not created by this function - the caller is responsible for creating it.
///
/// # Returns
///
/// A `PathBuf` representing the unique temporary directory path.
///
/// # Example
///
/// ```
/// let temp_dir = gdock::utils::get_unique_tempdir();
/// std::fs::create_dir_all(&temp_dir).expect("Failed to create temp dir");
/// // ... use the temp directory ...
/// std::fs::remove_dir_all(&temp_dir).ok();
/// ```
pub fn get_unique_tempdir() -> PathBuf {
    use std::env;
    use std::time::{SystemTime, UNIX_EPOCH};

    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .subsec_nanos();

    env::temp_dir().join(format!("gdock_{}_{}", std::process::id(), nanos))
}

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
    fn test_get_unique_tempdir() {
        let temp_dir = get_unique_tempdir();

        // Should be under system temp directory
        assert!(temp_dir.starts_with(std::env::temp_dir()));

        // Should contain "gdock_" and process ID
        let path_str = temp_dir.to_str().unwrap();
        assert!(path_str.contains("gdock_"));
        assert!(path_str.contains(&std::process::id().to_string()));
    }

    #[test]
    fn test_get_unique_tempdir_unique() {
        let temp_dir1 = get_unique_tempdir();
        let temp_dir2 = get_unique_tempdir();

        // Two consecutive calls should generate different paths
        assert_ne!(temp_dir1, temp_dir2);
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
