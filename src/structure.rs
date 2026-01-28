use crate::constants;
use crate::toppar;

use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct Model(pub Vec<Molecule>);

impl Default for Model {
    fn default() -> Self {
        Self::new()
    }
}

impl Model {
    pub fn new() -> Model {
        Model(Vec::new())
    }
}

#[derive(Debug, Clone)]
pub struct Molecule(pub Vec<Atom>);

#[derive(Debug, Clone)]
pub struct Atom {
    pub serial: i32,
    pub name: String,
    pub altloc: char,
    pub resname: String,
    pub chainid: char,
    pub resseq: i16,
    pub icode: char,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f32,
    pub tempfactor: f32,
    pub element: String,
    pub charge: f64,
    pub vdw_radius: f64,
    // pub sigma: f64,
    pub epsilon: f64,
    pub rmin2: f64,
    pub eps_1_4: f64,
    pub rmin2_1_4: f64,
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

impl Molecule {
    pub fn new() -> Molecule {
        Molecule(Vec::new())
    }

    // Function to compute the center of mass of the molecule
    pub fn center_of_mass(&self) -> (f64, f64, f64) {
        let num_atoms = self.0.len() as f64;
        let sum_x = self.0.iter().map(|atom| atom.x).sum::<f64>();
        let sum_y = self.0.iter().map(|atom| atom.y).sum::<f64>();
        let sum_z = self.0.iter().map(|atom| atom.z).sum::<f64>();
        (sum_x / num_atoms, sum_y / num_atoms, sum_z / num_atoms)
    }

    pub fn translate(&mut self, dx: f64, dy: f64, dz: f64) {
        for atom in &mut self.0 {
            atom.x += dx;
            atom.y += dy;
            atom.z += dz;
        }
    }

    // Function to rotate the molecule around its own center using Euler angles
    pub fn rotate(mut self, alpha: f64, beta: f64, gamma: f64) -> Self {
        // Translate the molecule to the center of mass
        let (center_x, center_y, center_z) = self.center_of_mass();
        for atom in &mut self.0 {
            atom.x -= center_x;
            atom.y -= center_y;
            atom.z -= center_z;
        }

        // Apply rotation using Euler angles
        for atom in &mut self.0 {
            let x = atom.x;
            let y = atom.y;
            let z = atom.z;

            // Apply rotation transformation
            let new_x = x * alpha.cos() * beta.cos()
                + y * (-alpha.sin() * gamma.cos() + alpha.cos() * beta.sin() * gamma.sin())
                + z * (alpha.sin() * gamma.sin() + alpha.cos() * beta.sin() * gamma.cos());
            let new_y = x * alpha.sin() * beta.cos()
                + y * (alpha.cos() * gamma.cos() + alpha.sin() * beta.sin() * gamma.sin())
                + z * (-alpha.cos() * gamma.sin() + alpha.sin() * beta.sin() * gamma.cos());
            let new_z =
                -x * beta.sin() + y * beta.cos() * gamma.sin() + z * beta.cos() * gamma.cos();

            atom.x = new_x;
            atom.y = new_y;
            atom.z = new_z;
        }

        // Move the molecule back to its original position
        for atom in &mut self.0 {
            atom.x += center_x;
            atom.y += center_y;
            atom.z += center_z;
        }
        self
    }

    pub fn displace(
        mut self,
        displacement_x: f64,
        displacement_y: f64,
        displacement_z: f64,
    ) -> Self {
        for atom in &mut self.0 {
            atom.x += displacement_x;
            atom.y += displacement_y;
            atom.z += displacement_z;
        }
        self
    }
}

impl Atom {
    fn new() -> Atom {
        Atom {
            serial: 0,
            name: String::new(),
            altloc: ' ',
            resname: String::new(),
            chainid: ' ',
            resseq: 0,
            icode: ' ',
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 0.0,
            tempfactor: 0.0,
            element: String::new(),
            charge: 0.0,
            vdw_radius: 0.0,
            epsilon: 0.0,
            rmin2: 0.0,
            eps_1_4: 0.0,
            rmin2_1_4: 0.0,
        }
    }
    // Function to convert Atom to PDB string format
    fn to_pdb_string(&self) -> String {
        format!(
            "ATOM  {:>5} {:^4}{:>1}{:3} {:>1}{:>4}{:>1}   {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}\n",
            self.serial,
            self.name,
            self.altloc,
            self.resname,
            self.chainid,
            self.resseq,
            self.icode,
            self.x,
            self.y,
            self.z,
            self.occupancy,
            self.tempfactor,
            self.element
        )
    }

    pub fn distance_to(&self, other: &Atom) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

/// Parse an ATOM line from a PDB file using fixed-column format
fn process_atom_line(line: &str) -> Option<Atom> {
    // ATOM lines must be at least 54 characters (up to Z coordinate)
    if !line.starts_with("ATOM") || line.len() < 54 {
        return None;
    }

    // Helper to safely extract and parse a substring
    let get_str = |start: usize, end: usize| -> &str {
        if line.len() >= end {
            &line[start..end]
        } else {
            ""
        }
    };

    // Parse using PDB fixed-column positions (1-indexed in spec, 0-indexed in code)
    let mut atom = Atom::new();

    // Columns 7-11: Atom serial number
    atom.serial = get_str(6, 11).trim().parse().ok()?;

    // Columns 13-16: Atom name
    atom.name = get_str(12, 16).trim().to_string();

    // Column 17: Alternate location indicator
    atom.altloc = get_str(16, 17).chars().next().unwrap_or(' ');

    // Columns 18-20: Residue name
    atom.resname = get_str(17, 20).trim().to_string();

    // Column 22: Chain identifier
    atom.chainid = get_str(21, 22).chars().next().unwrap_or(' ');

    // Columns 23-26: Residue sequence number
    atom.resseq = get_str(22, 26).trim().parse().ok()?;

    // Column 27: Code for insertions
    atom.icode = get_str(26, 27).chars().next().unwrap_or(' ');

    // Columns 31-38: X coordinate
    atom.x = get_str(30, 38).trim().parse().ok()?;

    // Columns 39-46: Y coordinate
    atom.y = get_str(38, 46).trim().parse().ok()?;

    // Columns 47-54: Z coordinate
    atom.z = get_str(46, 54).trim().parse().ok()?;

    // Columns 55-60: Occupancy (optional)
    atom.occupancy = get_str(54, 60).trim().parse().unwrap_or(1.0);

    // Columns 61-66: Temperature factor (optional)
    atom.tempfactor = get_str(60, 66).trim().parse().unwrap_or(0.0);

    // Columns 77-78: Element symbol (optional)
    atom.element = get_str(76, 78).trim().to_string();

    // If element is not specified, infer from atom name
    if atom.element.is_empty() {
        atom.element = atom.name.chars().next().unwrap_or(' ').to_string();
    }

    // Set VDW radius based on element
    atom.vdw_radius = match atom.element.trim() {
        "H" => constants::HYDROGEN_RADIUS,
        "C" => constants::CARBON_RADIUS,
        "N" => constants::NITROGEN_RADIUS,
        "O" => constants::OXYGEN_RADIUS,
        _ => 1.0,
    };

    // Get force field parameters from topology
    let atom_type = toppar::get_atom(atom.resname.as_str(), atom.name.as_str());

    if let Some(v) = atom_type {
        atom.epsilon = toppar::get_epsilon(v).unwrap_or(0.0);
        atom.rmin2 = toppar::get_rmin2(v).unwrap_or(0.0);
        atom.eps_1_4 = toppar::get_eps_1_4(v).unwrap_or(0.0);
        atom.rmin2_1_4 = toppar::get_rmin2_1_4(v).unwrap_or(0.0);
        atom.charge = toppar::get_charge(v).unwrap_or(0.0);
    }

    Some(atom)
}

// Reads a PDB file and returns a `Model` struct.
pub fn read_pdb(pdb_file: &str) -> Model {
    let mut model = Model::new();
    let mut molecule = Molecule::new();
    let mut has_model_records = false;

    let file = File::open(pdb_file).expect("Cannot open file");

    for line in BufReader::new(file).lines().map_while(Result::ok) {
        if line.starts_with("MODEL") {
            has_model_records = true;
            // Start a new molecule for this MODEL
            molecule = Molecule::new();
        } else if line.starts_with("ENDMDL") {
            // End of current model - push it to the model list
            model.0.push(molecule.clone());
            molecule = Molecule::new();
        } else if let Some(atom) = process_atom_line(&line) {
            molecule.0.push(atom);
        }
    }

    // If no MODEL records were found, push the single molecule
    if !has_model_records {
        model.0.push(molecule);
    }

    model
}

// Output a Molecule in PDB format
pub fn write_pdb(molecule: &Molecule, output_file: &str) {
    // Write the output string to the file
    let mut pdb_string = String::new();
    for atom in &molecule.0 {
        pdb_string.push_str(&atom.to_pdb_string());
    }
    fs::write(output_file, pdb_string).expect("Unable to write file");
}

pub fn distance(atom1: &Atom, atom2: &Atom) -> f64 {
    let dx = atom1.x - atom2.x;
    let dy = atom1.y - atom2.y;
    let dz = atom1.z - atom2.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

pub fn filter_by_resseq_vec(molecule: &Molecule, resseq_vec: &HashSet<i16>) -> Molecule {
    let mut filtered_molecule = Molecule::new();
    for atom in &molecule.0 {
        if resseq_vec.contains(&atom.resseq) {
            filtered_molecule.0.push(atom.clone());
        }
    }
    filtered_molecule
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn create_test_atom(x: f64, y: f64, z: f64) -> Atom {
        Atom {
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
    fn test_read_pdb() {
        let molecule = read_pdb(&"data/2oob.pdb".to_string());
        assert!(!molecule.0.is_empty())
    }

    #[test]
    fn test_molecule_new() {
        let mol = Molecule::new();
        assert_eq!(mol.0.len(), 0);
    }

    #[test]
    fn test_molecule_default() {
        let mol = Molecule::default();
        assert_eq!(mol.0.len(), 0);
    }

    #[test]
    fn test_center_of_mass_single_atom() {
        let mut mol = Molecule::new();
        mol.0.push(create_test_atom(1.0, 2.0, 3.0));

        let (cx, cy, cz) = mol.center_of_mass();
        assert_eq!(cx, 1.0);
        assert_eq!(cy, 2.0);
        assert_eq!(cz, 3.0);
    }

    #[test]
    fn test_center_of_mass_multiple_atoms() {
        let mut mol = Molecule::new();
        mol.0.push(create_test_atom(0.0, 0.0, 0.0));
        mol.0.push(create_test_atom(2.0, 0.0, 0.0));
        mol.0.push(create_test_atom(0.0, 2.0, 0.0));
        mol.0.push(create_test_atom(0.0, 0.0, 2.0));

        let (cx, cy, cz) = mol.center_of_mass();
        assert!((cx - 0.5).abs() < 1e-10);
        assert!((cy - 0.5).abs() < 1e-10);
        assert!((cz - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_translate() {
        let mut mol = Molecule::new();
        mol.0.push(create_test_atom(1.0, 2.0, 3.0));
        mol.0.push(create_test_atom(4.0, 5.0, 6.0));

        mol.translate(1.0, -1.0, 2.0);

        assert_eq!(mol.0[0].x, 2.0);
        assert_eq!(mol.0[0].y, 1.0);
        assert_eq!(mol.0[0].z, 5.0);
        assert_eq!(mol.0[1].x, 5.0);
        assert_eq!(mol.0[1].y, 4.0);
        assert_eq!(mol.0[1].z, 8.0);
    }

    #[test]
    fn test_displace() {
        let mut mol = Molecule::new();
        mol.0.push(create_test_atom(1.0, 2.0, 3.0));
        mol.0.push(create_test_atom(4.0, 5.0, 6.0));

        let mol = mol.displace(1.0, -1.0, 2.0);

        assert_eq!(mol.0[0].x, 2.0);
        assert_eq!(mol.0[0].y, 1.0);
        assert_eq!(mol.0[0].z, 5.0);
        assert_eq!(mol.0[1].x, 5.0);
        assert_eq!(mol.0[1].y, 4.0);
        assert_eq!(mol.0[1].z, 8.0);
    }

    #[test]
    fn test_rotate_zero_angles() {
        let mut mol = Molecule::new();
        mol.0.push(create_test_atom(1.0, 0.0, 0.0));

        let rotated = mol.rotate(0.0, 0.0, 0.0);

        // With zero rotation, position should be unchanged
        assert!((rotated.0[0].x - 1.0).abs() < 1e-10);
        assert!((rotated.0[0].y - 0.0).abs() < 1e-10);
        assert!((rotated.0[0].z - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_rotate_around_center() {
        let mut mol = Molecule::new();
        // Create a molecule with center at origin
        mol.0.push(create_test_atom(1.0, 0.0, 0.0));
        mol.0.push(create_test_atom(-1.0, 0.0, 0.0));

        // Rotate 90 degrees around z-axis
        let rotated = mol.rotate(0.0, 0.0, PI / 2.0);

        // After rotation, the center should still be at origin
        let (cx, cy, cz) = rotated.center_of_mass();
        assert!(cx.abs() < 1e-10);
        assert!(cy.abs() < 1e-10);
        assert!(cz.abs() < 1e-10);
    }

    #[test]
    fn test_atom_distance_to() {
        let atom1 = create_test_atom(0.0, 0.0, 0.0);
        let atom2 = create_test_atom(3.0, 4.0, 0.0);

        let dist = atom1.distance_to(&atom2);
        assert!((dist - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_distance_function() {
        let atom1 = create_test_atom(0.0, 0.0, 0.0);
        let atom2 = create_test_atom(1.0, 1.0, 1.0);

        let dist = distance(&atom1, &atom2);
        assert!((dist - 3.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_filter_by_resseq_vec() {
        let mut mol = Molecule::new();

        let mut atom1 = create_test_atom(0.0, 0.0, 0.0);
        atom1.resseq = 1;
        let mut atom2 = create_test_atom(1.0, 0.0, 0.0);
        atom2.resseq = 2;
        let mut atom3 = create_test_atom(2.0, 0.0, 0.0);
        atom3.resseq = 3;
        let mut atom4 = create_test_atom(3.0, 0.0, 0.0);
        atom4.resseq = 4;

        mol.0.push(atom1);
        mol.0.push(atom2);
        mol.0.push(atom3);
        mol.0.push(atom4);

        let mut filter_set = HashSet::new();
        filter_set.insert(2);
        filter_set.insert(4);

        let filtered = filter_by_resseq_vec(&mol, &filter_set);

        assert_eq!(filtered.0.len(), 2);
        assert_eq!(filtered.0[0].resseq, 2);
        assert_eq!(filtered.0[1].resseq, 4);
    }

    #[test]
    fn test_filter_by_resseq_vec_empty() {
        let mut mol = Molecule::new();
        mol.0.push(create_test_atom(0.0, 0.0, 0.0));

        let filter_set = HashSet::new();
        let filtered = filter_by_resseq_vec(&mol, &filter_set);

        assert_eq!(filtered.0.len(), 0);
    }

    #[test]
    fn test_read_pdb_single_model_no_model_record() {
        // Test with the actual single-model file in the data directory
        let pdb_path = "data/1h590A_0A.pdb".to_string();

        if !std::path::Path::new(&pdb_path).exists() {
            println!("Skipping test: {} not found", pdb_path);
            return;
        }

        let model = read_pdb(&pdb_path);

        // Should have exactly 1 molecule (no MODEL records)
        assert_eq!(model.0.len(), 1, "Single-model PDB should have 1 molecule");

        // The molecule should have many atoms
        assert!(
            !model.0[0].0.is_empty(),
            "Single molecule should have atoms"
        );

        println!("Single-model file has {} atoms", model.0[0].0.len());
    }

    #[test]
    fn test_read_pdb_multi_model_with_model_records() {
        // Test with the actual multi-model file in the data directory
        let pdb_path = "data/1h590B_0A_poses.pdb".to_string();

        if !std::path::Path::new(&pdb_path).exists() {
            println!("Skipping test: {} not found", pdb_path);
            return;
        }

        let model = read_pdb(&pdb_path);

        // Should have multiple molecules (has MODEL records)
        assert!(
            model.0.len() > 1,
            "Multi-model PDB should have more than 1 molecule, found {}",
            model.0.len()
        );

        // Each model should have atoms
        for (i, molecule) in model.0.iter().enumerate() {
            println!("Model {} has {} atoms", i + 1, molecule.0.len());
            assert!(
                !molecule.0.is_empty(),
                "Model {} should have atoms (total models: {})",
                i + 1,
                model.0.len()
            );
        }

        // All models should have the same number of atoms (same protein, different poses)
        let first_model_atoms = model.0[0].0.len();
        for (i, molecule) in model.0.iter().enumerate() {
            assert_eq!(
                molecule.0.len(),
                first_model_atoms,
                "Model {} has {} atoms, expected {}",
                i + 1,
                molecule.0.len(),
                first_model_atoms
            );
        }

        println!(
            "Multi-model file has {} models with {} atoms each",
            model.0.len(),
            first_model_atoms
        );
    }

    #[test]
    fn test_read_pdb_backward_compatibility() {
        // Verify that the test file used in other tests (data/2oob.pdb) works correctly
        let pdb_path = "data/2oob.pdb".to_string();

        if !std::path::Path::new(&pdb_path).exists() {
            println!("Skipping test: {} not found", pdb_path);
            return;
        }

        let model = read_pdb(&pdb_path);

        // This is a reference complex file - should be single model
        assert_eq!(model.0.len(), 1, "Reference complex should be single model");

        // Should have atoms
        assert!(
            !model.0[0].0.is_empty(),
            "Reference complex should have atoms"
        );

        println!("Reference complex file has {} atoms", model.0[0].0.len());
    }
}
