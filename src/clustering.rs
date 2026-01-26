//! FCC (Fraction of Common Contacts) clustering for protein structures.
//!
//! This module implements clustering based on contact similarity, adapted from
//! the HADDOCK FCC algorithm. Structures are grouped by the fraction of
//! inter-chain contacts they share.

use std::collections::{HashMap, HashSet};

use rayon::prelude::*;

use crate::structure::Molecule;

/// Configuration parameters for FCC clustering.
#[derive(Debug, Clone)]
pub struct ClusteringConfig {
    /// Distance cutoff for defining contacts (default: 5.0 Å)
    pub contact_distance: f64,
    /// Minimum FCC value to consider two structures as neighbors (default: 0.60)
    pub fcc_cutoff: f64,
    /// Strictness factor for asymmetric FCC comparison (default: 0.75)
    pub strictness: f64,
    /// Minimum cluster size to include in results (default: 4)
    pub min_cluster_size: usize,
}

impl Default for ClusteringConfig {
    fn default() -> Self {
        Self {
            contact_distance: 5.0,
            fcc_cutoff: 0.60,
            strictness: 0.75,
            min_cluster_size: 4,
        }
    }
}

/// Result of clustering a single group of structures.
#[derive(Debug, Clone)]
pub struct ClusterResult {
    /// Index of the cluster center in the original input
    pub center_idx: usize,
    /// Indices of all cluster members (including center)
    pub members: Vec<usize>,
    /// Number of members in this cluster
    pub size: usize,
}

/// Internal representation of a structure for clustering.
struct Element {
    neighbors: HashSet<usize>,
}

/// Type alias for residue coordinate map (chain, resseq) -> list of (x, y, z) coordinates
type ResidueCoordMap = HashMap<(char, i16), Vec<(f64, f64, f64)>>;

/// Calculate inter-chain contacts for a structure.
///
/// A contact is defined as any pair of residues from different chains
/// where at least one pair of heavy atoms is within the contact distance.
///
/// Returns a set of contact strings in format "chainA resA chainB resB".
pub fn calculate_contacts(molecule: &Molecule, contact_distance: f64) -> HashSet<String> {
    // Group atoms by residue (chain + resseq)
    let mut residues: ResidueCoordMap = HashMap::new();

    for atom in &molecule.0 {
        // Skip hydrogen atoms
        if atom.element.trim() == "H" {
            continue;
        }

        residues
            .entry((atom.chainid, atom.resseq))
            .or_default()
            .push((atom.x, atom.y, atom.z));
    }

    let mut contacts = HashSet::new();
    let residue_keys: Vec<_> = residues.keys().collect();

    // Compare all residue pairs from different chains
    for i in 0..residue_keys.len() {
        for j in (i + 1)..residue_keys.len() {
            let (chain_a, res_a) = residue_keys[i];
            let (chain_b, res_b) = residue_keys[j];

            // Only consider inter-chain contacts
            if chain_a == chain_b {
                continue;
            }

            let atoms_a = &residues[residue_keys[i]];
            let atoms_b = &residues[residue_keys[j]];

            // Check if any atom pair is within contact distance
            'outer: for (xa, ya, za) in atoms_a {
                for (xb, yb, zb) in atoms_b {
                    let dist = ((xa - xb).powi(2) + (ya - yb).powi(2) + (za - zb).powi(2)).sqrt();
                    if dist <= contact_distance {
                        // Store contact in canonical order (smaller chain first)
                        let contact = if chain_a < chain_b {
                            format!("{} {} {} {}", chain_a, res_a, chain_b, res_b)
                        } else {
                            format!("{} {} {} {}", chain_b, res_b, chain_a, res_a)
                        };
                        contacts.insert(contact);
                        break 'outer;
                    }
                }
            }
        }
    }

    contacts
}

/// Calculate the Fraction of Common Contacts between two contact sets.
///
/// Returns (fcc, fcc_v) where:
/// - fcc = |X ∩ Y| / |X|
/// - fcc_v = |X ∩ Y| / |Y|
fn calculate_fcc(x: &HashSet<String>, y: &HashSet<String>) -> (f64, f64) {
    if x.is_empty() || y.is_empty() {
        return (0.0, 0.0);
    }

    let common = x.intersection(y).count() as f64;
    let fcc = common / x.len() as f64;
    let fcc_v = common / y.len() as f64;

    (fcc, fcc_v)
}

/// Calculate pairwise FCC values for all structure pairs (parallel).
fn calculate_pairwise_fcc(contact_sets: &[HashSet<String>]) -> Vec<(usize, usize, f64, f64)> {
    let n = contact_sets.len();

    // Generate all (i, j) pairs and compute FCC in parallel
    (0..n)
        .into_par_iter()
        .flat_map(|i| {
            (0..n).into_par_iter().map(move |j| {
                let (fcc, fcc_v) = calculate_fcc(&contact_sets[i], &contact_sets[j]);
                (i, j, fcc, fcc_v)
            })
        })
        .collect()
}

/// Build the similarity graph from pairwise FCC values.
fn create_elements(
    pairwise_fcc: Vec<(usize, usize, f64, f64)>,
    n_structures: usize,
    config: &ClusteringConfig,
) -> HashMap<usize, Element> {
    let mut elements: HashMap<usize, Element> = HashMap::new();

    // Initialize all elements
    for i in 0..n_structures {
        elements.insert(
            i,
            Element {
                neighbors: HashSet::new(),
            },
        );
    }

    // Add neighbor relationships based on FCC threshold
    for (i, j, fcc, fcc_v) in pairwise_fcc {
        // Bidirectional threshold check
        if fcc >= config.fcc_cutoff && fcc_v >= config.fcc_cutoff * config.strictness {
            elements.get_mut(&i).unwrap().neighbors.insert(j);
        }
        if fcc_v >= config.fcc_cutoff && fcc >= config.fcc_cutoff * config.strictness {
            elements.get_mut(&j).unwrap().neighbors.insert(i);
        }
    }

    elements
}

/// Perform greedy clustering based on the similarity graph.
fn cluster_elements(mut elements: HashMap<usize, Element>) -> Vec<ClusterResult> {
    let mut used: HashSet<usize> = HashSet::new();
    let mut clusters = Vec::new();

    loop {
        // Get clusterable elements (not yet used)
        let clusterable: Vec<usize> = elements
            .keys()
            .filter(|k| !used.contains(k))
            .copied()
            .collect();

        if clusterable.is_empty() {
            break;
        }

        // Find the element with the most neighbors (greedy selection)
        // Tie-break by index (lower index preferred for determinism)
        let center = clusterable
            .iter()
            .max_by(|&&a, &&b| {
                let count_a = elements[&a]
                    .neighbors
                    .iter()
                    .filter(|n| !used.contains(n))
                    .count();
                let count_b = elements[&b]
                    .neighbors
                    .iter()
                    .filter(|n| !used.contains(n))
                    .count();
                count_a.cmp(&count_b).then_with(|| b.cmp(&a))
            })
            .copied()
            .unwrap();

        // Get neighbors that aren't used yet
        let neighbors: Vec<usize> = elements[&center]
            .neighbors
            .iter()
            .filter(|n| !used.contains(n))
            .copied()
            .collect();

        // Build cluster members (center + neighbors)
        let mut members: Vec<usize> = vec![center];
        for neighbor in &neighbors {
            if *neighbor != center {
                members.push(*neighbor);
                used.insert(*neighbor);
            }
        }
        used.insert(center);

        // Sort members for consistent output
        members.sort();

        clusters.push(ClusterResult {
            center_idx: center,
            members: members.clone(),
            size: members.len(),
        });

        // Remove center from elements
        elements.remove(&center);
    }

    // Sort clusters by size (descending)
    clusters.sort_by(|a, b| b.size.cmp(&a.size));

    clusters
}

/// Cluster structures based on contact similarity (FCC).
///
/// # Arguments
/// * `structures` - Slice of molecules (each should be a receptor+ligand complex)
/// * `config` - Clustering configuration parameters
///
/// # Returns
/// Vector of `ClusterResult` sorted by cluster size (largest first).
/// Only clusters meeting `min_cluster_size` are returned.
pub fn cluster_structures(
    structures: &[Molecule],
    config: &ClusteringConfig,
) -> Vec<ClusterResult> {
    if structures.is_empty() {
        return Vec::new();
    }

    // Calculate contacts for each structure (parallel)
    let contact_sets: Vec<HashSet<String>> = structures
        .par_iter()
        .map(|mol| calculate_contacts(mol, config.contact_distance))
        .collect();

    // Calculate pairwise FCC
    let pairwise_fcc = calculate_pairwise_fcc(&contact_sets);

    // Build similarity graph
    let elements = create_elements(pairwise_fcc, structures.len(), config);

    // Perform clustering
    let clusters = cluster_elements(elements);

    // Filter by minimum cluster size
    clusters
        .into_iter()
        .filter(|c| c.size >= config.min_cluster_size)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clustering_config_default() {
        let config = ClusteringConfig::default();
        assert!((config.contact_distance - 5.0).abs() < 0.001);
        assert!((config.fcc_cutoff - 0.60).abs() < 0.001);
        assert!((config.strictness - 0.75).abs() < 0.001);
        assert_eq!(config.min_cluster_size, 4);
    }

    #[test]
    fn test_calculate_fcc_identical() {
        let mut set_a = HashSet::new();
        set_a.insert("A 1 B 2".to_string());
        set_a.insert("A 3 B 4".to_string());

        let (fcc, fcc_v) = calculate_fcc(&set_a, &set_a);
        assert!((fcc - 1.0).abs() < 0.001);
        assert!((fcc_v - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_calculate_fcc_disjoint() {
        let mut set_a = HashSet::new();
        set_a.insert("A 1 B 2".to_string());

        let mut set_b = HashSet::new();
        set_b.insert("A 5 B 6".to_string());

        let (fcc, fcc_v) = calculate_fcc(&set_a, &set_b);
        assert!((fcc - 0.0).abs() < 0.001);
        assert!((fcc_v - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_calculate_fcc_partial_overlap() {
        let mut set_a = HashSet::new();
        set_a.insert("A 1 B 2".to_string());
        set_a.insert("A 3 B 4".to_string());

        let mut set_b = HashSet::new();
        set_b.insert("A 1 B 2".to_string());
        set_b.insert("A 5 B 6".to_string());

        let (fcc, fcc_v) = calculate_fcc(&set_a, &set_b);
        assert!((fcc - 0.5).abs() < 0.001);
        assert!((fcc_v - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_calculate_fcc_empty() {
        let set_a: HashSet<String> = HashSet::new();
        let set_b: HashSet<String> = HashSet::new();

        let (fcc, fcc_v) = calculate_fcc(&set_a, &set_b);
        assert!((fcc - 0.0).abs() < 0.001);
        assert!((fcc_v - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_cluster_structures_empty() {
        let structures: Vec<Molecule> = Vec::new();
        let config = ClusteringConfig::default();
        let clusters = cluster_structures(&structures, &config);
        assert!(clusters.is_empty());
    }

    #[test]
    fn test_calculate_contacts_with_real_data() {
        use crate::commands::run::combine_molecules;
        use crate::structure::read_pdb;

        let receptor_model = read_pdb(&"data/A.pdb".to_string());
        let ligand_model = read_pdb(&"data/B.pdb".to_string());

        let receptor = &receptor_model.0[0];
        let ligand = &ligand_model.0[0];

        let complex = combine_molecules(receptor, ligand);
        let contacts = calculate_contacts(&complex, 5.0);

        // Should have inter-chain contacts (A-B or B-A)
        assert!(!contacts.is_empty(), "Should detect inter-chain contacts");

        // All contacts should involve two different chains
        for contact in &contacts {
            let parts: Vec<&str> = contact.split_whitespace().collect();
            assert_eq!(
                parts.len(),
                4,
                "Contact format should be 'chain res chain res'"
            );
            assert_ne!(parts[0], parts[2], "Contacts should be inter-chain");
        }
    }
}
