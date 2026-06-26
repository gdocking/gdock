use crate::structure::{read_pdb, Atom, Molecule};
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

/// Generate restraints from interface contacts and print to stdout.
///
/// Uses HADDOCK-style cross-product: finds all receptor residues with any contact
/// to the ligand (active rec) and all ligand residues with any contact to the
/// receptor (active lig), then emits every active_rec × active_lig pair.
/// Each anchor gets the full active ligand set as its OR group, matching the
/// ambiguity semantics of HADDOCK AIR restraints.
pub fn generate_restraints(receptor_file: String, ligand_file: String, cutoff: f64) {
    let receptor_model = read_pdb(&receptor_file);
    let ligand_model = read_pdb(&ligand_file);

    let receptor = &receptor_model.0[0];
    let ligand = &ligand_model.0[0];

    let interface_pairs = find_interface_pairs(receptor, ligand, cutoff);

    let mut active_rec: Vec<i16> = interface_pairs.iter().map(|(r, _)| *r).collect();
    active_rec.sort();
    active_rec.dedup();

    let mut active_lig: Vec<i16> = interface_pairs.iter().map(|(_, l)| *l).collect();
    active_lig.sort();
    active_lig.dedup();

    let mut pairs: Vec<String> = Vec::new();
    for r in &active_rec {
        for l in &active_lig {
            pairs.push(format!("{}:{}", r, l));
        }
    }

    println!("{}", pairs.join(","));
}

/// Minimum heavy-atom distance from a pre-collected receptor atom slice to a ligand residue.
fn min_dist_to_res(rec_atoms: &[&Atom], ligand: &Molecule, lig_res: i16) -> f64 {
    ligand
        .0
        .iter()
        .filter(|a| a.resseq == lig_res && !a.name.starts_with('H'))
        .flat_map(|la| {
            rec_atoms.iter().map(move |ra| {
                let dx = ra.x - la.x;
                let dy = ra.y - la.y;
                let dz = ra.z - la.z;
                (dx * dx + dy * dy + dz * dz).sqrt()
            })
        })
        .fold(f64::INFINITY, f64::min)
}

/// Like `find_interface_pairs` but returns only the single closest ligand partner
/// per receptor residue (minimum heavy-atom distance).
///
/// Use this to generate unambiguous 1:1 restraints where the GA must bring a
/// specific receptor residue close to a specific ligand residue, rather than
/// satisfying any member of a large OR group.
pub fn find_unambiguous_pairs(
    receptor: &Molecule,
    ligand: &Molecule,
    cutoff: f64,
) -> Vec<(i16, i16)> {
    let all_pairs = find_interface_pairs(receptor, ligand, cutoff);

    let mut by_rec: HashMap<i16, Vec<i16>> = HashMap::new();
    for (r, l) in &all_pairs {
        by_rec.entry(*r).or_default().push(*l);
    }

    let mut result: Vec<(i16, i16)> = by_rec
        .iter()
        .filter_map(|(rec_res, lig_candidates)| {
            let rec_atoms: Vec<&Atom> = receptor
                .0
                .iter()
                .filter(|a| a.resseq == *rec_res && !a.name.starts_with('H'))
                .collect();
            lig_candidates
                .iter()
                .map(|&l| (l, min_dist_to_res(&rec_atoms, ligand, l)))
                .min_by(|(_, d1), (_, d2)| d1.partial_cmp(d2).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(closest, _)| (*rec_res, closest))
        })
        .collect();

    result.sort();
    result
}

/// Generate unambiguous restraints: one closest-partner pair per receptor residue.
///
/// Unlike `generate_restraints` (which emits a cross-product of all active residues),
/// this selects only the single ligand residue that is closest (minimum heavy-atom
/// distance) to each receptor interface residue. The resulting pairs are tight
/// 1:1 restraints suitable for AmbiguousRestraint with a single-element OR group.
pub fn generate_restraints_unambig(receptor_file: String, ligand_file: String, cutoff: f64) {
    let receptor_model = read_pdb(&receptor_file);
    let ligand_model = read_pdb(&ligand_file);

    let receptor = &receptor_model.0[0];
    let ligand = &ligand_model.0[0];

    let unambig_pairs = find_unambiguous_pairs(receptor, ligand, cutoff);

    let pairs: Vec<String> = unambig_pairs
        .iter()
        .map(|(r, l)| format!("{}:{}", r, l))
        .collect();

    println!("{}", pairs.join(","));
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

/// Parse a HADDOCK AIR `.tbl` file into active-active `(receptor_resseq, ligand_resseq)` pairs.
///
/// Chain identities are detected automatically: the segid that first appears as an
/// anchor is treated as the receptor; the other segid is the ligand.  Only residues
/// that appear as anchors are considered active; passive OR-group members (never
/// anchors) are silently dropped.
pub fn tbl_to_pairs(tbl_file: &str) -> Vec<(i32, i32)> {
    let content = std::fs::read_to_string(tbl_file).unwrap_or_else(|e| {
        eprintln!("Error: cannot read TBL file '{}': {}", tbl_file, e);
        std::process::exit(1);
    });

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

    let sel_re = Regex::new(
        r"(?i)\(\s*(?:resid?\s+(\d+)\s+and\s+segid\s+(\w+)|segid\s+(\w+)\s+and\s+resid?\s+(\d+))\s*\)",
    )
    .expect("invalid regex");

    let assign_re = Regex::new(r"(?i)\bassign\b").expect("invalid regex");
    let blocks: Vec<&str> = assign_re.split(&stripped).skip(1).collect();

    let parsed: Vec<Vec<(String, i32)>> = blocks
        .iter()
        .map(|b| parse_selections(b, &sel_re))
        .collect();

    // Auto-detect chains: first segid seen as an anchor = receptor; the other = ligand.
    let mut rec = String::new();
    let mut lig = String::new();
    'outer: for sels in &parsed {
        if let Some((anchor_seg, _)) = sels.first() {
            if rec.is_empty() {
                rec = anchor_seg.clone();
            }
            for (seg, _) in &sels[1..] {
                if seg != &rec && lig.is_empty() {
                    lig = seg.clone();
                    break 'outer;
                }
            }
        }
    }

    if rec.is_empty() || lig.is_empty() {
        eprintln!(
            "Error: could not detect two distinct chain segids in '{}'",
            tbl_file
        );
        std::process::exit(1);
    }

    // First pass — every anchor is active.
    let mut active: HashSet<(String, i32)> = HashSet::new();
    for sels in &parsed {
        if let Some(anchor) = sels.first() {
            active.insert(anchor.clone());
        }
    }

    // Second pass — collect pairs where both anchor and partner are active.
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

    let mut pairs: Vec<(i32, i32)> = groups
        .into_iter()
        .flat_map(|(rec_resi, lig_set)| lig_set.into_iter().map(move |l| (rec_resi, l)))
        .collect();
    pairs.sort();

    if pairs.is_empty() {
        eprintln!(
            "Warning: no active-active pairs found for chains '{}'/'{}'",
            rec, lig
        );
    }

    pairs
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

    #[test]
    fn test_find_unambiguous_pairs_one_per_receptor_residue() {
        let receptor_model = read_pdb("data/2oob_A.pdb");
        let ligand_model = read_pdb("data/2oob_B.pdb");

        let receptor = &receptor_model.0[0];
        let ligand = &ligand_model.0[0];

        let unambig = find_unambiguous_pairs(receptor, ligand, 5.0);

        // Each receptor residue must appear exactly once
        let mut rec_residues: Vec<i16> = unambig.iter().map(|(r, _)| *r).collect();
        let original_len = rec_residues.len();
        rec_residues.dedup();
        assert_eq!(
            rec_residues.len(),
            original_len,
            "Each receptor residue should appear at most once"
        );

        // Result must be a subset of find_interface_pairs (every pair must be a valid contact)
        let all_pairs = find_interface_pairs(receptor, ligand, 5.0);
        for (r, l) in &unambig {
            assert!(
                all_pairs.contains(&(*r, *l)),
                "Unambig pair ({r},{l}) must be a valid interface contact"
            );
        }
    }

    #[test]
    fn test_find_unambiguous_pairs_fewer_than_ambig() {
        let receptor_model = read_pdb("data/2oob_A.pdb");
        let ligand_model = read_pdb("data/2oob_B.pdb");

        let receptor = &receptor_model.0[0];
        let ligand = &ligand_model.0[0];

        let ambig = find_interface_pairs(receptor, ligand, 5.0);
        let unambig = find_unambiguous_pairs(receptor, ligand, 5.0);

        // Unambig yields at most one pair per receptor residue, so never more than ambig
        assert!(
            unambig.len() <= ambig.len(),
            "Unambig pairs ({}) should be <= interface pairs ({})",
            unambig.len(),
            ambig.len()
        );
        assert!(!unambig.is_empty(), "Should find at least one unambig pair");
    }

    #[test]
    fn test_find_unambiguous_pairs_selects_closest() {
        use crate::structure::{Atom, Molecule};

        fn make_atom(resseq: i16, x: f64, y: f64, z: f64) -> Atom {
            Atom {
                serial: 1,
                name: "CA".to_string(),
                altloc: ' ',
                resname: "ALA".to_string(),
                chainid: 'A',
                resseq,
                icode: ' ',
                x,
                y,
                z,
                occupancy: 1.0,
                tempfactor: 0.0,
                element: "C".to_string(),
                charge: 0.0,
                vdw_radius: 1.7,
                epsilon: -0.1,
                rmin2: 2.0,
                eps_1_4: -0.1,
                rmin2_1_4: 1.9,
            }
        }

        // Receptor residue 1 at origin; two ligand residues: 10 at 3Å, 11 at 8Å.
        // Both within cutoff=10. Unambig must pick residue 10 (closer).
        let mut receptor = Molecule::new();
        receptor.0.push(make_atom(1, 0.0, 0.0, 0.0));

        let mut ligand = Molecule::new();
        ligand.0.push(make_atom(10, 3.0, 0.0, 0.0));
        ligand.0.push(make_atom(11, 8.0, 0.0, 0.0));

        let pairs = find_unambiguous_pairs(&receptor, &ligand, 10.0);
        assert_eq!(pairs.len(), 1);
        assert_eq!(
            pairs[0],
            (1, 10),
            "Should select residue 10 (3Å) over 11 (8Å)"
        );
    }
}
