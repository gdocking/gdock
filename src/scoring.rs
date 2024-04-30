use crate::structure;
use regex::Regex;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

// Ugly functions related to scoring
/// Reads a complex from a PDB file and returns two separate molecules.
///
/// # Arguments
///
/// * `pdbf` - A string slice that represents the path to the PDB file.
///
/// # Returns
///
/// A tuple containing two `structure::Molecule` objects representing the two separate molecules
/// extracted from the complex.
///
/// # Panics
///
/// This function will panic if it fails to remove the temporary files created during the process.
pub fn read_complex(pdbf: &str) -> (structure::Molecule, structure::Molecule) {
    // The input is a complex with multiple chains
    let molecules = split_complex(pdbf);

    let mol_a = structure::read_pdb(&molecules[0]);
    let mol_b = structure::read_pdb(&molecules[1]);

    // Delete the files
    std::fs::remove_file(&molecules[0]).unwrap();
    std::fs::remove_file(&molecules[1]).unwrap();

    (mol_a, mol_b)
}

// Read a complex PDB and write its receptor and ligand to separate PDB files
/// Splits a PDB file into separate files based on the chain identifier.
///
/// This function takes a PDB file path as input and reads the atoms from the file.
/// It then populates a hashmap with atoms of each chain, where the chain identifier
/// is extracted from the atom line using a regular expression.
/// Finally, it creates separate PDB files for each chain, containing the atoms
/// corresponding to that chain, and returns a vector of the file names.
///
/// # Arguments
///
/// * `pdb_file` - A string slice representing the path to the PDB file.
///
/// # Returns
///
/// A vector of strings representing the file names of the split PDB files.
pub fn split_complex(pdb_file: &str) -> Vec<String> {
    // Populate the hashmap with atoms of a given chain
    let mut atom_map: HashMap<String, Vec<String>> = HashMap::new();

    // Open the complex_file and read the atoms
    // let file = File::open(complex_file).expect("Cannot open file");
    let file = File::open(pdb_file).expect("Cannot open file");
    let pdb_re = Regex::new(r"ATOM\s{2}((?:\s|\d){5})\s((?:\s|.){4})((?:\s|\w){1})((?:\s|\w){3})\s((?:\s|\w){1})((?:\s|\w){4})((?:\s|\w){1})\s{3}((?:\s|.){8})((?:\s|.){8})((?:\s|.){8})((?:\s|.){6})((?:\s|.){6})\s{10}((?:\s|\w){2})((?:\s|\w){2})").unwrap();
    for line in BufReader::new(file).lines().flatten() {
        if let Some(cap) = pdb_re.captures(&line) {
            let chain = cap[5].trim().to_string();

            // Add it to the hashmap if it doesn't exist
            if !atom_map.contains_key(&chain) {
                atom_map.insert(chain.clone(), Vec::new());
            }

            // Add the line to the hashmap with a newline
            let mut atom = line.clone();
            atom.push('\n');
            atom_map.get_mut(&chain).unwrap().push(atom);
        }
    }

    let mut result = vec![];

    for (chain, atoms) in atom_map {
        let fname = format!("{}.pdb", chain);
        let mut file = File::create(&fname).expect("Cannot create file");
        let mut pdb_string = String::new();

        for atom in atoms {
            pdb_string.push_str(&atom);
        }

        file.write_all(pdb_string.as_bytes())
            .expect("Cannot write to file");

        // Add the filename to the result vector
        result.push(fname);
    }
    result
}
