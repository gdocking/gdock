use crate::structure;
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

    // NOTE: Use only the first models
    let mol_a = structure::read_pdb(&molecules[0]).0[0].clone();
    let mol_b = structure::read_pdb(&molecules[1]).0[0].clone();

    // Delete the files and temp directory
    std::fs::remove_file(&molecules[0]).unwrap();
    std::fs::remove_file(&molecules[1]).unwrap();
    if let Some(parent) = std::path::Path::new(&molecules[0]).parent() {
        let _ = std::fs::remove_dir(parent);
    }

    (mol_a, mol_b)
}

// Read a complex PDB and write its receptor and ligand to separate PDB files
/// Splits a PDB file into separate files based on the chain identifier.
///
/// This function takes a PDB file path as input and reads the atoms from the file.
/// It then populates a hashmap with atoms of each chain, where the chain identifier
/// is extracted from the atom line using a regular expression.
/// Finally, it creates separate PDB files for each chain in a temporary directory,
/// containing the atoms corresponding to that chain, and returns a vector of the
/// absolute file paths.
///
/// # Arguments
///
/// * `pdb_file` - A string slice representing the path to the PDB file.
///
/// # Returns
///
/// A vector of strings representing the absolute file paths of the split PDB files
/// in the temporary directory.
pub fn split_complex(pdb_file: &str) -> Vec<String> {
    // Create a unique temporary directory for this operation
    let temp_dir = crate::utils::get_unique_tempdir();
    std::fs::create_dir_all(&temp_dir).expect("Cannot create temp directory");

    // Populate the hashmap with atoms of a given chain
    let mut atom_map: HashMap<String, Vec<String>> = HashMap::new();

    // Open the complex_file and read the atoms
    let file = File::open(pdb_file).expect("Cannot open file");
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        // Parse ATOM lines by column position (PDB is fixed-column format)
        // Chain ID is at column 22 (0-indexed: 21)
        if line.starts_with("ATOM") && line.len() >= 22 {
            let chain = line
                .chars()
                .nth(21)
                .unwrap_or(' ')
                .to_string()
                .trim()
                .to_string();
            let chain = if chain.is_empty() {
                " ".to_string()
            } else {
                chain
            };

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

    // Sort chains alphabetically to ensure deterministic order
    let mut chains: Vec<_> = atom_map.keys().cloned().collect();
    chains.sort();

    let mut result = vec![];

    for chain in chains {
        let atoms = &atom_map[&chain];
        let file_path = temp_dir.join(format!("{}.pdb", chain));
        let mut file = File::create(&file_path).expect("Cannot create file");
        let mut pdb_string = String::new();

        for atom in atoms {
            pdb_string.push_str(atom);
        }

        file.write_all(pdb_string.as_bytes())
            .expect("Cannot write to file");

        // Add the absolute filename to the result vector
        result.push(file_path.to_str().unwrap().to_string());
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;

    #[test]
    fn test_split_complex_creates_two_chains() {
        let result = split_complex("data/2oob.pdb");

        // Should have 2 chains (A and B)
        assert_eq!(result.len(), 2);

        // Verify both files exist
        assert!(Path::new(&result[0]).exists());
        assert!(Path::new(&result[1]).exists());

        // Verify files are in temp directory
        for path in &result {
            assert!(path.contains("gdock_"));
        }

        // Clean up
        if let Some(parent) = Path::new(&result[0]).parent() {
            let _ = fs::remove_dir_all(parent);
        }
    }

    #[test]
    fn test_split_complex_chain_names() {
        let result = split_complex("data/2oob.pdb");

        // Check that files are named A.pdb and B.pdb
        assert!(result[0].ends_with("A.pdb"));
        assert!(result[1].ends_with("B.pdb"));

        // Clean up
        if let Some(parent) = Path::new(&result[0]).parent() {
            let _ = fs::remove_dir_all(parent);
        }
    }

    #[test]
    fn test_split_complex_files_have_content() {
        let result = split_complex("data/2oob.pdb");

        // Read and verify both files have content
        let content_a = fs::read_to_string(&result[0]).expect("Failed to read chain A file");
        let content_b = fs::read_to_string(&result[1]).expect("Failed to read chain B file");

        assert!(!content_a.is_empty());
        assert!(!content_b.is_empty());

        // Verify they contain ATOM lines
        assert!(content_a.contains("ATOM"));
        assert!(content_b.contains("ATOM"));

        // Clean up
        if let Some(parent) = Path::new(&result[0]).parent() {
            let _ = fs::remove_dir_all(parent);
        }
    }

    #[test]
    fn test_read_complex_returns_two_molecules() {
        let (mol_a, mol_b) = read_complex("data/2oob.pdb");

        // Both molecules should have atoms
        assert!(!mol_a.0.is_empty());
        assert!(!mol_b.0.is_empty());

        // Chain A should have atoms with chain ID 'A'
        assert_eq!(mol_a.0[0].chainid, 'A');

        // Chain B should have atoms with chain ID 'B'
        assert_eq!(mol_b.0[0].chainid, 'B');
    }

    #[test]
    fn test_read_complex_cleans_up_temp_files() {
        // Call read_complex which should create and then delete temp files
        let (mol_a, mol_b) = read_complex("data/2oob.pdb");

        // Verify molecules were created
        assert!(!mol_a.0.is_empty());
        assert!(!mol_b.0.is_empty());

        // We can't easily verify cleanup without modifying the function,
        // but we can at least verify the function completed successfully
    }

    #[test]
    fn test_read_complex_preserves_coordinates() {
        let (mol_a, mol_b) = read_complex("data/2oob.pdb");

        // Verify coordinates are valid numbers (not NaN or Inf)
        for atom in &mol_a.0 {
            assert!(atom.x.is_finite());
            assert!(atom.y.is_finite());
            assert!(atom.z.is_finite());
        }

        for atom in &mol_b.0 {
            assert!(atom.x.is_finite());
            assert!(atom.y.is_finite());
            assert!(atom.z.is_finite());
        }
    }

    #[test]
    #[should_panic(expected = "Cannot open file")]
    fn test_split_complex_nonexistent_file() {
        split_complex("nonexistent_file_xyz123.pdb");
    }

    #[test]
    fn test_split_complex_deterministic_order() {
        // Run twice and verify we get the same order
        let result1 = split_complex("data/2oob.pdb");
        let result2 = split_complex("data/2oob.pdb");

        // Extract just the filenames (not full paths)
        let name1_a = Path::new(&result1[0]).file_name().unwrap();
        let name1_b = Path::new(&result1[1]).file_name().unwrap();
        let name2_a = Path::new(&result2[0]).file_name().unwrap();
        let name2_b = Path::new(&result2[1]).file_name().unwrap();

        // Should be in same order (alphabetically sorted)
        assert_eq!(name1_a, name2_a);
        assert_eq!(name1_b, name2_b);

        // Clean up
        if let Some(parent) = Path::new(&result1[0]).parent() {
            let _ = fs::remove_dir_all(parent);
        }
        if let Some(parent) = Path::new(&result2[0]).parent() {
            let _ = fs::remove_dir_all(parent);
        }
    }
}
