use crate::constants;
use crate::toppar;
use regex::Regex;

use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct Molecule(pub Vec<Atom>);

#[derive(Debug, Clone)]
pub struct Atom {
    pub serial: i32,
    pub name: String,
    altloc: char,
    resname: String,
    pub chainid: char,
    pub resseq: i16,
    icode: char,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    occupancy: f32,
    tempfactor: f32,
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

    // // Add the restraints to the molecule
    // pub fn add_restraints(&mut self, restraints: Vec<i16>) {
    //     self.1 = restraints;
    // }
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

/// Reads a PDB file and returns a `Molecule` struct.
///
/// # Arguments
///
/// * `pdb_file` - A reference to a `String` containing the path to the PDB file.
///
/// # Returns
///
/// A `Molecule` struct containing the atoms parsed from the PDB file.
pub fn read_pdb(pdb_file: &String) -> Molecule {
    let mut molecule = Molecule::new();

    // Read the PDB file
    let file = File::open(pdb_file).expect("Cannot open file");
    let pdb_re = Regex::new(r"ATOM\s{2}((?:\s|\d){5})\s((?:\s|.){4})((?:\s|\w){1})((?:\s|\w){3})\s((?:\s|\w){1})((?:\s|\w){4})((?:\s|\w){1})\s{3}((?:\s|.){8})((?:\s|.){8})((?:\s|.){8})((?:\s|.){6})((?:\s|.){6})\s{10}((?:\s|\w){2})((?:\s|\w){2})").unwrap();
    for line in BufReader::new(file).lines().flatten() {
        if let Some(cap) = pdb_re.captures(&line) {
            let mut atom = Atom::new();
            atom.serial = cap[1].trim().parse().unwrap();
            atom.name = cap[2].trim().to_string();
            atom.altloc = cap[3].parse().unwrap();
            atom.resname = cap[4].trim().to_string();
            atom.chainid = cap[5].trim().parse().unwrap();
            atom.resseq = cap[6].trim().parse().unwrap();
            atom.icode = cap[7].parse().unwrap();
            atom.x = cap[8].trim().parse().unwrap();
            atom.y = cap[9].trim().parse().unwrap();
            atom.z = cap[10].trim().parse().unwrap();
            atom.occupancy = cap[11].trim().parse().unwrap();
            atom.tempfactor = cap[12].trim().parse().unwrap();
            atom.element = cap[13].trim().to_string();
            // atom.charge = cap[14].trim().to_string();
            atom.vdw_radius = match atom.element.as_str() {
                "H" => constants::HYDROGEN_RADIUS,
                "C" => constants::CARBON_RADIUS,
                "N" => constants::NITROGEN_RADIUS,
                "O" => constants::OXYGEN_RADIUS,
                _ => 1.0,
            };

            let atom_type = toppar::get_atom(atom.resname.as_str(), atom.name.as_str());

            match atom_type {
                Some(v) => {
                    atom.epsilon = toppar::get_epsilon(v).unwrap();
                    atom.rmin2 = toppar::get_rmin2(v).unwrap();
                    atom.eps_1_4 = toppar::get_eps_1_4(v).unwrap();
                    atom.rmin2_1_4 = toppar::get_rmin2_1_4(v).unwrap();
                    atom.charge = toppar::get_charge(v).unwrap();
                }
                None => continue,
            }

            // let charge = toppar::get_charge(atom_type);

            molecule.0.push(atom);
        }
    }

    // // Read the restraint file and add the restraints to the molecule
    // let file = File::open(restraint_file).expect("Cannot open file");
    // let restraint_re = Regex::new(r"^(\d+)").unwrap();
    // for line in BufReader::new(file).lines().flatten() {
    //     if let Some(cap) = restraint_re.captures(&line) {
    //         let mut restraint = Restraint::new();
    //         restraint.resnum = cap[1].trim().parse().unwrap();
    //         molecule.1.push(restraint);
    //     }
    // }

    molecule
}

// Output a Molecule in PDB format
#[allow(dead_code)]
pub fn write_pdb(molecule: &Molecule, output_file: &String) {
    // let data = "Some data!";
    // fs::write("/tmp/foo", data).expect("Unable to write file");
    // }

    // let mut file = File::create(output_file).expect("Cannot create file");

    // Write the output string to the file
    let mut pdb_string = String::new();
    for atom in &molecule.0 {
        pdb_string.push_str(&atom.to_pdb_string());
        // let s = format!(
        //     "ATOM  {:>5} {:<4}{:>1}{:<3} {:>1}{:>4}{:>1}   {:>8.3} {:>8.3} {:>8.3}{:>6.2}{:>6.2}          {:>2}{:>2}\n",
        //     atom.serial, atom.name, atom.altloc, atom.resname, atom.chainud, atom.resseq, atom.icode, atom.x, atom.y, atom.z, atom.occupancy, atom.tempfactor, atom.element, atom.charge
        // );
        // file.write_all(pdb_string.as_bytes())
        //     .expect("Cannot write to file");
    }
    fs::write(output_file, pdb_string).expect("Unable to write file");
}

// Change the chain of a Molecule
#[allow(dead_code)]
fn change_chain(molecule: &mut Molecule, chain: char) {
    for atom in &mut molecule.0 {
        atom.chainid = chain;
    }
}

#[allow(dead_code)]
pub fn distance(atom1: &Atom, atom2: &Atom) -> f64 {
    let dx = atom1.x - atom2.x;
    let dy = atom1.y - atom2.y;
    let dz = atom1.z - atom2.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// fn generate_random_rotation() -> (f64, f64, f64) {
//     let mut rng = thread_rng();
//     let rotation_x = rng.gen_range(0.0..=2.0 * PI);
//     let rotation_y = rng.gen_range(0.0..=2.0 * PI);
//     let rotation_z = rng.gen_range(0.0..=2.0 * PI);
//     (rotation_x, rotation_y, rotation_z)
// }

// fn generate_random_displacement(max_displacement: f64) -> (f64, f64, f64) {
//     let mut rng = thread_rng();
//     let displacement_x = rng.gen_range(-max_displacement..=max_displacement);
//     let displacement_y = rng.gen_range(-max_displacement..=max_displacement);
//     let displacement_z = rng.gen_range(-max_displacement..=max_displacement);
//     (displacement_x, displacement_y, displacement_z)
// }

pub fn filter_by_resseq_vec(molecule: &Molecule, resseq_vec: &HashSet<i16>) -> Molecule {
    let mut filtered_molecule = Molecule::new();
    for atom in &molecule.0 {
        if resseq_vec.contains(&atom.resseq) {
            filtered_molecule.0.push(atom.clone());
        }
    }
    filtered_molecule
}
