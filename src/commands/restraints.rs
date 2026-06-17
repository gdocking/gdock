use crate::structure::{read_pdb, Molecule};
use regex::Regex;
use std::collections::{HashMap, HashSet};

/// Find interface residue pairs between receptor and ligand within a distance cutoff.
/// Returns pairs of (receptor_resseq, ligand_resseq) that are in contact.
pub fn find_interface_pairs(
    receptor: &Molecule,
    ligand: &Molecule,
    cutoff: f64,
) -> Vec<(i16, i16)> {
    let mut interface_pairs: Vec<(i16, i16)> = Vec::new();

    // Get unique residues from receptor and ligand
    let rec_residues: HashSet<i16> = receptor.0.iter().map(|a| a.resseq).collect();
    let lig_residues: HashSet<i16> = ligand.0.iter().map(|a| a.resseq).collect();

    for rec_res in &rec_residues {
        for lig_res in &lig_residues {
            // Get heavy atoms for this residue pair
            let rec_atoms: Vec<_> = receptor
                .0
                .iter()
                .filter(|a| a.resseq == *rec_res && !a.name.starts_with('H'))
                .collect();
            let lig_atoms: Vec<_> = ligand
                .0
                .iter()
                .filter(|a| a.resseq == *lig_res && !a.name.starts_with('H'))
                .collect();

            // Check if any pair of atoms is within cutoff
            'outer: for rec_atom in &rec_atoms {
                for lig_atom in &lig_atoms {
                    let dx = rec_atom.x - lig_atom.x;
                    let dy = rec_atom.y - lig_atom.y;
                    let dz = rec_atom.z - lig_atom.z;
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt();

                    if dist <= cutoff {
                        interface_pairs.push((*rec_res, *lig_res));
                        break 'outer;
                    }
                }
            }
        }
    }

    // Sort for consistent output
    interface_pairs.sort();
    interface_pairs
}

/// Generate restraints from interface contacts and print to stdout
pub fn generate_restraints(receptor_file: String, ligand_file: String, cutoff: f64) {
    let receptor_model = read_pdb(&receptor_file);
    let ligand_model = read_pdb(&ligand_file);

    let receptor = &receptor_model.0[0];
    let ligand = &ligand_model.0[0];

    let interface_pairs = find_interface_pairs(receptor, ligand, cutoff);

    // Output in gdock restraints format
    let restraints_str: Vec<String> = interface_pairs
        .iter()
        .map(|(r, l)| format!("{}:{}", r, l))
        .collect();

    println!("{}", restraints_str.join(","));
}

/// Extract `(segid, resseq)` tuples from one assign block's text.
///
/// Handles both `(resid N and segid X)` and `(segid X and resid N)` forms,
/// and accepts the `resi` abbreviation used by some HADDOCK versions.
fn parse_selections(block: &str, sel_re: &Regex) -> Vec<(String, i32)> {
    sel_re
        .captures_iter(block)
        .filter_map(|cap| {
            if let (Some(resi), Some(segid)) = (cap.get(1), cap.get(2)) {
                let n = resi.as_str().parse::<i32>().ok()?;
                Some((segid.as_str().to_uppercase(), n))
            } else if let (Some(segid), Some(resi)) = (cap.get(3), cap.get(4)) {
                let n = resi.as_str().parse::<i32>().ok()?;
                Some((segid.as_str().to_uppercase(), n))
            } else {
                None
            }
        })
        .collect()
}

/// Convert a HADDOCK AIR `.tbl` file to gdock active-active restraint pairs.
///
/// Only residues that appear as anchors (left-hand side of `assign` blocks) are
/// considered active.  A pair `receptor_resi:ligand_resi` is emitted only when
/// both the anchor and its OR-group partner are active.  Passive OR-group members
/// (never anchors) are silently dropped.
///
/// Output is written to stdout as a comma-separated string in `rec:lig` format,
/// ready to be passed to `gdock run --restraints`.
pub fn convert_tbl(tbl_file: &str, receptor_chain: &str, ligand_chain: &str) {
    let content = std::fs::read_to_string(tbl_file).unwrap_or_else(|e| {
        eprintln!("Error: cannot read TBL file '{}': {}", tbl_file, e);
        std::process::exit(1);
    });

    // Strip inline comments — everything after `!` on each line.
    let stripped: String = content
        .lines()
        .map(|line| {
            if let Some(pos) = line.find('!') {
                &line[..pos]
            } else {
                line
            }
        })
        .collect::<Vec<_>>()
        .join("\n");

    // Match (resid N and segid X) or (segid X and resid N); `resi` also accepted.
    let sel_re = Regex::new(
        r"(?i)\(\s*(?:resid?\s+(\d+)\s+and\s+segid\s+(\w+)|segid\s+(\w+)\s+and\s+resid?\s+(\d+))\s*\)",
    )
    .expect("invalid regex");

    // Split on the `assign` keyword to get one chunk per restraint.
    let assign_re = Regex::new(r"(?i)\bassign\b").expect("invalid regex");
    let blocks: Vec<&str> = assign_re.split(&stripped).skip(1).collect();

    let rec = receptor_chain.to_uppercase();
    let lig = ligand_chain.to_uppercase();

    // First pass — every anchor is active.
    let mut active: HashSet<(String, i32)> = HashSet::new();
    let parsed: Vec<Vec<(String, i32)>> = blocks
        .iter()
        .map(|b| parse_selections(b, &sel_re))
        .collect();
    for sels in &parsed {
        if let Some(anchor) = sels.first() {
            active.insert(anchor.clone());
        }
    }

    // Second pass — emit pairs where both anchor and partner are active.
    let mut groups: HashMap<i32, HashSet<i32>> = HashMap::new();
    for sels in &parsed {
        if sels.len() < 2 {
            continue;
        }
        let (anchor_seg, anchor_resi) = &sels[0];
        for (partner_seg, partner_resi) in &sels[1..] {
            if !active.contains(&(partner_seg.clone(), *partner_resi)) {
                continue;
            }
            if anchor_seg == &rec && partner_seg == &lig {
                groups
                    .entry(*anchor_resi)
                    .or_default()
                    .insert(*partner_resi);
            } else if anchor_seg == &lig && partner_seg == &rec {
                groups
                    .entry(*partner_resi)
                    .or_default()
                    .insert(*anchor_resi);
            }
        }
    }

    // Flatten to sorted pairs and print.
    let mut pairs: Vec<(i32, i32)> = groups
        .into_iter()
        .flat_map(|(rec_resi, lig_set)| lig_set.into_iter().map(move |l| (rec_resi, l)))
        .collect();
    pairs.sort();

    if pairs.is_empty() {
        eprintln!(
            "Warning: no active-active pairs found for chains '{}'/'{}'",
            receptor_chain, ligand_chain
        );
        return;
    }

    let output: Vec<String> = pairs.iter().map(|(r, l)| format!("{}:{}", r, l)).collect();
    println!("{}", output.join(","));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_interface_pairs_with_test_data() {
        let receptor_model = read_pdb("data/2oob_A.pdb");
        let ligand_model = read_pdb("data/2oob_B.pdb");

        let receptor = &receptor_model.0[0];
        let ligand = &ligand_model.0[0];

        let pairs = find_interface_pairs(receptor, ligand, 5.0);

        // Should find some interface pairs
        assert!(!pairs.is_empty(), "Should find interface pairs");

        // All pairs should be sorted
        let mut sorted = pairs.clone();
        sorted.sort();
        assert_eq!(pairs, sorted, "Pairs should be sorted");
    }

    #[test]
    fn test_find_interface_pairs_cutoff_effect() {
        let receptor_model = read_pdb("data/2oob_A.pdb");
        let ligand_model = read_pdb("data/2oob_B.pdb");

        let receptor = &receptor_model.0[0];
        let ligand = &ligand_model.0[0];

        let pairs_5a = find_interface_pairs(receptor, ligand, 5.0);
        let pairs_10a = find_interface_pairs(receptor, ligand, 10.0);

        // Larger cutoff should find more or equal pairs
        assert!(
            pairs_10a.len() >= pairs_5a.len(),
            "10Å cutoff should find >= pairs than 5Å"
        );
    }
}
