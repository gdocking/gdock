use crate::restraints;
use crate::structure;
use crate::structure::Atom;

// Calculate the ratio of clashes
//   i.e. the number of atom pairs that are closer than the sum of their van der Waals radii
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

/// Calculates the Lennard-Jones potential energy between two particles.
///
/// The Lennard-Jones potential energy is a pairwise interaction potential used to model
/// the intermolecular forces between particles. It is commonly used in molecular dynamics simulations.
fn lj_potential(atom1: &Atom, atom2: &Atom, distance: f64) -> f64 {
    let epsilon_ij = (atom1.epsilon * atom2.epsilon).sqrt();
    let rmin_ij = atom1.rmin2 + atom2.rmin2;
    let ratio_pow = (rmin_ij / distance).powi(6);
    epsilon_ij * (ratio_pow.powi(2) - 2.0 * ratio_pow)
}

/// Calculates the van der Waals energy between two molecules.
///
/// # Arguments
///
/// * `molecule1` - The first molecule.
/// * `molecule2` - The second molecule.
///
/// # Returns
///
/// The van der Waals energy between the two molecules.
pub fn vdw_energy(receptor: &structure::Molecule, ligand: &structure::Molecule) -> f64 {
    let mut energy = 0.0;
    let cutoff = 12.0; // VDW cutoff in Angstroms

    for atom1 in &receptor.0 {
        for atom2 in &ligand.0 {
            // Check if atoms are not from the same molecule
            if atom1.serial != atom2.serial {
                let dist = structure::distance(atom1, atom2);
                if dist > 0.0 && dist < cutoff {
                    energy += lj_potential(atom1, atom2, dist);
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
    let cutoff = 8.0; // Cutoff for burial calculation in Angstroms

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
            if dist < cutoff {
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
            if dist < cutoff {
                burial_count += 1;
            }
        }

        let burial = (burial_count as f64 / 10.0).min(1.0);
        energy -= asp2 * burial * atom2.vdw_radius * atom2.vdw_radius * 4.0 * std::f64::consts::PI;
    }

    energy
}

pub fn bsa_energy(_molecule1: &structure::Molecule, _molecule2: &structure::Molecule) -> f64 {
    todo!();
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
    let distance_threshold = 8.0; // Distance at which restraint is satisfied (matches restraint definition)
    let force_constant = 10.0; // kcal/mol per Å² - controls penalty strength

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

            // Harmonic penalty for violated restraints
            // Only applies when CA-CA distance > threshold
            if dist > distance_threshold {
                let violation = dist - distance_threshold;
                energy += force_constant * violation * violation;
            }
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
    let cutoff = 15.0; // Electrostatics cutoff in Angstroms
    let min_dist = 1.0; // Minimum distance to avoid singularities

    for atom1 in &receptor.0 {
        for atom2 in &ligand.0 {
            let dist = structure::distance(atom1, atom2);
            if dist > min_dist && dist < cutoff {
                let charge1 = atom1.charge;
                let charge2 = atom2.charge;
                // Distance-dependent dielectric: ε = r
                energy += k * (charge1 * charge2) / (dist * dist);
            }
        }
    }

    energy
}
