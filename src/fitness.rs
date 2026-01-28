use crate::constants::{
    AIR_FORCE_CONSTANT, DESOLV_CUTOFF, ELEC_CUTOFF, ELEC_MIN_DISTANCE, SOFTCORE_ALPHA, VDW_CUTOFF,
};
use crate::restraints;
use crate::structure;
use crate::structure::Atom;

/// Calculates the ratio of clashing atom pairs between two molecules.
///
/// A clash is defined as an atom pair where the distance is less than
/// the sum of their van der Waals radii.
///
/// # Arguments
///
/// * `molecule1` - First molecule
/// * `molecule2` - Second molecule
///
/// # Returns
///
/// The ratio of clashing pairs to total pairs (0.0 to 1.0).
pub fn calc_clashes(molecule1: &structure::Molecule, molecule2: &structure::Molecule) -> f64 {
    let mut clashes = 0;
    let mut total_atoms = 0;

    for atom1 in &molecule1.0 {
        for atom2 in &molecule2.0 {
            total_atoms += 1;

            let dist = structure::distance(atom1, atom2);
            let vdw_dist = atom1.vdw_radius + atom2.vdw_radius;
            if dist < vdw_dist {
                clashes += 1;
            }
        }
    }
    clashes as f64 / total_atoms as f64
}

/// Soft-core Lennard-Jones potential for docking.
///
/// Uses a modified LJ potential that remains finite at r→0, allowing
/// gentle overlap without catastrophic energies. This enables the GA
/// to explore near-native conformations that may have minor clashes
/// before refinement.
///
/// The soft-core potential transitions smoothly:
/// - At large distances (r > rmin): standard LJ behavior
/// - At close distances (r < 0.8*rmin): soft repulsion with finite maximum
///
/// # Arguments
/// * `atom1`, `atom2` - The two atoms
/// * `distance` - Distance between atoms in Angstroms
/// * `alpha` - Softness parameter (default 0.5, higher = softer)
///
/// # Returns
/// Soft-core VDW energy in kcal/mol
fn softcore_lj_potential(atom1: &Atom, atom2: &Atom, distance: f64, alpha: f64) -> f64 {
    let epsilon_ij = (atom1.epsilon * atom2.epsilon).sqrt();
    let rmin_ij = atom1.rmin2 + atom2.rmin2;

    // Soft-core modification: r_eff^6 = r^6 + alpha*rmin^6
    // This prevents singularity at r→0
    let rmin6 = rmin_ij.powi(6);
    let r6 = distance.powi(6);
    let r_eff6 = r6 + alpha * rmin6;

    // Standard LJ form but with effective distance
    let ratio = rmin6 / r_eff6;
    epsilon_ij * (ratio.powi(2) - 2.0 * ratio)
}

/// Calculates the van der Waals energy between two molecules.
///
/// Uses soft-core VDW potential to prevent catastrophic energies from clashes
/// while still penalizing severe overlaps. This allows the GA to explore
/// near-native conformations that may have minor clashes.
///
/// # Arguments
///
/// * `receptor` - The receptor molecule.
/// * `ligand` - The ligand molecule.
///
/// # Returns
///
/// The van der Waals energy between the two molecules in kcal/mol.
pub fn vdw_energy(receptor: &structure::Molecule, ligand: &structure::Molecule) -> f64 {
    let mut energy = 0.0;

    for atom1 in &receptor.0 {
        for atom2 in &ligand.0 {
            // Check if atoms are not from the same molecule
            if atom1.serial != atom2.serial {
                let dist = structure::distance(atom1, atom2);
                if dist > 0.0 && dist < VDW_CUTOFF {
                    // Use soft-core potential to prevent catastrophic energies
                    let vdw = softcore_lj_potential(atom1, atom2, dist, SOFTCORE_ALPHA);
                    // Still apply capping for extreme cases
                    let capped_vdw = vdw.clamp(-50.0, 500.0);
                    energy += capped_vdw;
                }
            }
        }
    }
    energy
}

/// Calculates the desolvation energy between two molecules.
///
/// Uses empirical atomic solvation parameters (ASP) to estimate the
/// desolvation penalty when atoms become buried upon complex formation.
/// This is crucial for protein-protein docking as it compensates for
/// the loss of favorable water interactions.
///
/// # Arguments
///
/// * `receptor` - The receptor molecule.
/// * `ligand` - The ligand molecule.
///
/// # Returns
///
/// The desolvation energy in kcal/mol (positive = unfavorable).
pub fn desolv_energy(receptor: &structure::Molecule, ligand: &structure::Molecule) -> f64 {
    let mut energy = 0.0;

    // Empirical atomic solvation parameters (kcal/mol/Å²)
    // Based on OPLS/CHARMM atom types - simplified version
    // Positive values = hydrophobic (burial is favorable)
    // Negative values = hydrophilic (burial is unfavorable)
    let get_asp = |element: &str| -> f64 {
        match element.trim() {
            "C" => 0.012,  // Carbon: hydrophobic, burial favorable
            "N" => -0.160, // Nitrogen: polar, burial unfavorable
            "O" => -0.160, // Oxygen: polar, burial unfavorable
            "S" => 0.012,  // Sulfur: somewhat hydrophobic
            "H" => 0.0,    // Hydrogen: usually ignored
            _ => 0.0,
        }
    };

    // Calculate desolvation for receptor atoms buried by ligand
    for atom1 in &receptor.0 {
        let asp1 = get_asp(&atom1.element);
        if asp1.abs() < 0.001 {
            continue; // Skip atoms with no ASP
        }

        let mut burial_count = 0;
        for atom2 in &ligand.0 {
            let dist = structure::distance(atom1, atom2);
            if dist < DESOLV_CUTOFF {
                burial_count += 1;
            }
        }

        // Burial factor: 0 = not buried, 1 = fully buried
        // Use sigmoid-like function for smooth burial
        let burial = (burial_count as f64 / 10.0).min(1.0);

        // Desolvation = ASP * burial * surface_area
        // Surface area of sphere = 4πr²
        energy -= asp1 * burial * atom1.vdw_radius * atom1.vdw_radius * 4.0 * std::f64::consts::PI;
    }

    // Calculate desolvation for ligand atoms buried by receptor
    for atom2 in &ligand.0 {
        let asp2 = get_asp(&atom2.element);
        if asp2.abs() < 0.001 {
            continue;
        }

        let mut burial_count = 0;
        for atom1 in &receptor.0 {
            let dist = structure::distance(atom1, atom2);
            if dist < DESOLV_CUTOFF {
                burial_count += 1;
            }
        }

        let burial = (burial_count as f64 / 10.0).min(1.0);
        energy -= asp2 * burial * atom2.vdw_radius * atom2.vdw_radius * 4.0 * std::f64::consts::PI;
    }

    energy
}

/// Calculates AIR (Ambiguous Interaction Restraints) energy.
///
/// AIR energy uses restraints as soft distance constraints to guide docking
/// toward native-like conformations. This is the core of information-driven
/// docking approaches like HADDOCK.
///
/// For each restraint that is NOT satisfied (distance > threshold), we add
/// a penalty that grows with distance. This creates a smooth energy landscape
/// that guides the search toward satisfying restraints.
///
/// # Arguments
///
/// * `restraints` - Vector of restraints from native structure
/// * `receptor` - The receptor molecule
/// * `ligand` - The ligand molecule
///
/// # Returns
///
/// The AIR energy penalty in kcal/mol (0 = all satisfied, positive = violations).
pub fn air_energy(
    restraints: &[crate::restraints::Restraint],
    receptor: &structure::Molecule,
    ligand: &structure::Molecule,
) -> f64 {
    let mut energy = 0.0;

    // HADDOCK-style flat-bottom restraints with lower/upper bounds
    // For CA-CA interface contacts (typical at interface: 3-7 Å)
    let lower_bound = 0.0; // No penalty for being close (contacts are good!)
    let upper_bound = 7.0; // Target maximum for interface CA-CA distance

    for restraint in restraints {
        // Find CA atoms for this restraint (consistent with is_satisfied)
        let ca_receptor = receptor
            .0
            .iter()
            .find(|a| a.resseq == restraint.0.resseq && a.name.trim() == "CA");

        let ca_ligand = ligand
            .0
            .iter()
            .find(|a| a.resseq == restraint.1.resseq && a.name.trim() == "CA");

        // Both CAs must exist
        if let (Some(ca1), Some(ca2)) = (ca_receptor, ca_ligand) {
            let dist = structure::distance(ca1, ca2);

            // Flat-bottom harmonic potential (HADDOCK-style):
            // No penalty if lower_bound <= dist <= upper_bound
            // Quadratic penalty outside this range
            let violation = if dist < lower_bound {
                lower_bound - dist // Too close (rarely happens for CA-CA)
            } else if dist > upper_bound {
                dist - upper_bound // Too far (main penalty)
            } else {
                0.0 // Within bounds - no penalty
            };

            energy += AIR_FORCE_CONSTANT * violation * violation;
        }
    }

    energy
}

pub fn satisfaction_ratio(
    restraints: &[restraints::Restraint],
    receptor: &structure::Molecule,
    ligand: &structure::Molecule,
) -> f64 {
    let satisfied = restraints
        .iter()
        .filter(|&x| x.is_satisfied(receptor, ligand))
        .count();

    let total = restraints.len();
    if total > 0 {
        satisfied as f64 / total as f64
    } else {
        0.0
    }
}

/// Calculates the electrostatic energy between two molecules.
///
/// Uses distance-dependent dielectric (ε = r) which is common in protein docking
/// to dampen long-range electrostatic interactions.
///
/// Formula: E = k * (q1 * q2) / (ε * r) where ε = r, thus E = k * (q1 * q2) / r²
///
/// # Arguments
///
/// * `receptor` - The receptor molecule.
/// * `ligand` - The ligand molecule.
///
/// # Returns
///
/// The electrostatic energy in kcal/mol.
pub fn elec_energy(receptor: &structure::Molecule, ligand: &structure::Molecule) -> f64 {
    // Coulomb's constant in kcal⋅Å/(mol⋅e²) for molecular mechanics
    let k = 332.0636;
    let mut energy = 0.0;

    for atom1 in &receptor.0 {
        for atom2 in &ligand.0 {
            let dist = structure::distance(atom1, atom2);
            if dist > ELEC_MIN_DISTANCE && dist < ELEC_CUTOFF {
                let charge1 = atom1.charge;
                let charge2 = atom2.charge;
                // Distance-dependent dielectric: ε = r
                energy += k * (charge1 * charge2) / (dist * dist);
            }
        }
    }

    energy
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_atom(x: f64, y: f64, z: f64, element: &str, charge: f64) -> structure::Atom {
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
            element: element.to_string(),
            charge,
            vdw_radius: 1.7,
            epsilon: -0.1,
            rmin2: 2.0,
            eps_1_4: -0.1,
            rmin2_1_4: 1.9,
        }
    }

    #[test]
    fn test_calc_clashes_no_clashes() {
        let mut mol1 = structure::Molecule::new();
        mol1.0.push(create_test_atom(0.0, 0.0, 0.0, "C", 0.0));

        let mut mol2 = structure::Molecule::new();
        mol2.0.push(create_test_atom(10.0, 0.0, 0.0, "C", 0.0)); // Far away

        let ratio = calc_clashes(&mol1, &mol2);
        assert_eq!(ratio, 0.0);
    }

    #[test]
    fn test_calc_clashes_with_clashes() {
        let mut mol1 = structure::Molecule::new();
        mol1.0.push(create_test_atom(0.0, 0.0, 0.0, "C", 0.0));

        let mut mol2 = structure::Molecule::new();
        mol2.0.push(create_test_atom(1.0, 0.0, 0.0, "C", 0.0)); // Within VDW distance (< 3.4Å)

        let ratio = calc_clashes(&mol1, &mol2);
        assert!(ratio > 0.0, "Should detect clash");
        assert!(ratio <= 1.0, "Clash ratio should be <= 1.0");
    }

    #[test]
    fn test_vdw_energy_far_atoms() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, "C", 0.0));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom(20.0, 0.0, 0.0, "C", 0.0)); // Beyond cutoff (12Å)

        let energy = vdw_energy(&receptor, &ligand);
        assert_eq!(energy, 0.0, "Energy should be zero beyond cutoff");
    }

    #[test]
    fn test_vdw_energy_close_atoms() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, "C", 0.0));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom(5.0, 0.0, 0.0, "C", 0.0)); // Within cutoff

        let energy = vdw_energy(&receptor, &ligand);
        assert!(energy.is_finite(), "Energy should be finite");
    }

    #[test]
    fn test_vdw_energy_very_close_atoms() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, "C", 0.0));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom(1.0, 0.0, 0.0, "C", 0.0)); // Very close (clash)

        let energy = vdw_energy(&receptor, &ligand);
        // Soft-core potential should keep this finite
        assert!(energy.is_finite(), "Soft-core should prevent infinities");
        assert!(energy <= 500.0, "Energy should be capped at 500.0");
    }

    #[test]
    fn test_elec_energy_neutral_atoms() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, "C", 0.0)); // No charge

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom(5.0, 0.0, 0.0, "C", 0.0)); // No charge

        let energy = elec_energy(&receptor, &ligand);
        assert_eq!(energy, 0.0, "Energy should be zero for neutral atoms");
    }

    #[test]
    fn test_elec_energy_charged_atoms() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, "C", 1.0)); // +1 charge

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom(5.0, 0.0, 0.0, "C", -1.0)); // -1 charge

        let energy = elec_energy(&receptor, &ligand);
        assert!(
            energy < 0.0,
            "Opposite charges should give attractive (negative) energy"
        );
        assert!(energy.is_finite(), "Energy should be finite");
    }

    #[test]
    fn test_elec_energy_same_charge() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, "C", 1.0)); // +1 charge

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom(5.0, 0.0, 0.0, "C", 1.0)); // +1 charge

        let energy = elec_energy(&receptor, &ligand);
        assert!(
            energy > 0.0,
            "Same charges should give repulsive (positive) energy"
        );
    }

    #[test]
    fn test_desolv_energy_hydrophobic_burial() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, "C", 0.0)); // Hydrophobic carbon

        let mut ligand = structure::Molecule::new();
        // Add many atoms nearby to create burial
        for i in 0..15 {
            ligand.0.push(create_test_atom(
                3.0 * (i as f64).cos(),
                3.0 * (i as f64).sin(),
                0.0,
                "C",
                0.0,
            ));
        }

        let energy = desolv_energy(&receptor, &ligand);
        // Hydrophobic burial should be favorable (negative energy)
        assert!(
            energy < 0.0,
            "Burying hydrophobic carbon should be favorable"
        );
    }

    #[test]
    fn test_desolv_energy_hydrophilic_burial() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, "O", 0.0)); // Hydrophilic oxygen

        let mut ligand = structure::Molecule::new();
        // Add many atoms nearby to create burial
        for i in 0..15 {
            ligand.0.push(create_test_atom(
                3.0 * (i as f64).cos(),
                3.0 * (i as f64).sin(),
                0.0,
                "C",
                0.0,
            ));
        }

        let energy = desolv_energy(&receptor, &ligand);
        // Hydrophilic burial should be unfavorable (positive energy)
        assert!(
            energy > 0.0,
            "Burying hydrophilic oxygen should be unfavorable"
        );
    }

    #[test]
    fn test_air_energy_satisfied_restraints() {
        let mut receptor = structure::Molecule::new();
        let mut atom1 = create_test_atom(0.0, 0.0, 0.0, "C", 0.0);
        atom1.name = "CA".to_string();
        atom1.resseq = 1;
        receptor.0.push(atom1.clone());

        let mut ligand = structure::Molecule::new();
        let mut atom2 = create_test_atom(5.0, 0.0, 0.0, "C", 0.0); // Within 7Å
        atom2.name = "CA".to_string();
        atom2.resseq = 10;
        ligand.0.push(atom2.clone());

        let restraints = vec![restraints::Restraint(atom1, atom2)];

        let energy = air_energy(&restraints, &receptor, &ligand);
        assert_eq!(energy, 0.0, "Satisfied restraints should have zero penalty");
    }

    #[test]
    fn test_air_energy_violated_restraints() {
        let mut receptor = structure::Molecule::new();
        let mut atom1 = create_test_atom(0.0, 0.0, 0.0, "C", 0.0);
        atom1.name = "CA".to_string();
        atom1.resseq = 1;
        receptor.0.push(atom1.clone());

        let mut ligand = structure::Molecule::new();
        let mut atom2 = create_test_atom(15.0, 0.0, 0.0, "C", 0.0); // Beyond 7Å
        atom2.name = "CA".to_string();
        atom2.resseq = 10;
        ligand.0.push(atom2.clone());

        let restraints = vec![restraints::Restraint(atom1, atom2)];

        let energy = air_energy(&restraints, &receptor, &ligand);
        assert!(energy > 0.0, "Violated restraints should have penalty");
    }

    #[test]
    fn test_satisfaction_ratio_all_satisfied() {
        let mut receptor = structure::Molecule::new();
        let mut atom1 = create_test_atom(0.0, 0.0, 0.0, "C", 0.0);
        atom1.name = "CA".to_string();
        atom1.resseq = 1;
        receptor.0.push(atom1.clone());

        let mut ligand = structure::Molecule::new();
        let mut atom2 = create_test_atom(5.0, 0.0, 0.0, "C", 0.0);
        atom2.name = "CA".to_string();
        atom2.resseq = 10;
        ligand.0.push(atom2.clone());

        let restraints = vec![restraints::Restraint(atom1, atom2)];

        let ratio = satisfaction_ratio(&restraints, &receptor, &ligand);
        assert_eq!(ratio, 1.0, "All restraints satisfied should give ratio 1.0");
    }

    #[test]
    fn test_satisfaction_ratio_none_satisfied() {
        let mut receptor = structure::Molecule::new();
        let mut atom1 = create_test_atom(0.0, 0.0, 0.0, "C", 0.0);
        atom1.name = "CA".to_string();
        atom1.resseq = 1;
        receptor.0.push(atom1.clone());

        let mut ligand = structure::Molecule::new();
        let mut atom2 = create_test_atom(15.0, 0.0, 0.0, "C", 0.0);
        atom2.name = "CA".to_string();
        atom2.resseq = 10;
        ligand.0.push(atom2.clone());

        let restraints = vec![restraints::Restraint(atom1, atom2)];

        let ratio = satisfaction_ratio(&restraints, &receptor, &ligand);
        assert_eq!(ratio, 0.0, "No restraints satisfied should give ratio 0.0");
    }

    #[test]
    fn test_satisfaction_ratio_empty_restraints() {
        let receptor = structure::Molecule::new();
        let ligand = structure::Molecule::new();
        let restraints = vec![];

        let ratio = satisfaction_ratio(&restraints, &receptor, &ligand);
        assert_eq!(ratio, 0.0, "Empty restraints should give ratio 0.0");
    }

    #[test]
    fn test_satisfaction_ratio_partial() {
        let mut receptor = structure::Molecule::new();
        let mut atom1 = create_test_atom(0.0, 0.0, 0.0, "C", 0.0);
        atom1.name = "CA".to_string();
        atom1.resseq = 1;
        let mut atom3 = create_test_atom(0.0, 0.0, 0.0, "C", 0.0);
        atom3.name = "CA".to_string();
        atom3.resseq = 2;
        receptor.0.push(atom1.clone());
        receptor.0.push(atom3.clone());

        let mut ligand = structure::Molecule::new();
        let mut atom2 = create_test_atom(5.0, 0.0, 0.0, "C", 0.0); // Satisfied
        atom2.name = "CA".to_string();
        atom2.resseq = 10;
        let mut atom4 = create_test_atom(15.0, 0.0, 0.0, "C", 0.0); // Violated
        atom4.name = "CA".to_string();
        atom4.resseq = 11;
        ligand.0.push(atom2.clone());
        ligand.0.push(atom4.clone());

        let restraints = vec![
            restraints::Restraint(atom1, atom2),
            restraints::Restraint(atom3, atom4),
        ];

        let ratio = satisfaction_ratio(&restraints, &receptor, &ligand);
        assert_eq!(ratio, 0.5, "1 of 2 satisfied should give ratio 0.5");
    }
}
