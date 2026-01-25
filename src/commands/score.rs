use crate::constants::EnergyWeights;
use crate::evaluator;
use crate::fitness;
use crate::restraints;
use crate::scoring;
use crate::structure::read_pdb;

/// Score structures and print results to stdout in TSV format
pub fn score(
    receptor_file: String,
    ligand_file: String,
    restraint_pairs: Option<Vec<(i32, i32)>>,
    reference_file: Option<String>,
    weights: EnergyWeights,
) {
    let receptor_model = read_pdb(&receptor_file);
    let ligand_model = read_pdb(&ligand_file);

    let receptor = receptor_model.0[0].clone();

    // Print header line for parser-friendly output
    // clash_pct is always included for calibration purposes
    if restraint_pairs.is_some() && reference_file.is_some() {
        println!("model\tscore\tvdw\telec\tdesolv\tair\tclash_pct\tw_vdw\tw_elec\tw_desolv\tw_air\tdockq\tlrmsd\tirmsd\tfnat");
    } else if restraint_pairs.is_some() {
        println!("model\tscore\tvdw\telec\tdesolv\tair\tclash_pct\tw_vdw\tw_elec\tw_desolv\tw_air");
    } else if reference_file.is_some() {
        println!("model\tscore\tvdw\telec\tdesolv\tclash_pct\tw_vdw\tw_elec\tw_desolv\tw_air\tdockq\tlrmsd\tirmsd\tfnat");
    } else {
        println!("model\tscore\tvdw\telec\tdesolv\tclash_pct\tw_vdw\tw_elec\tw_desolv\tw_air");
    }

    for (model_idx, ligand) in ligand_model.0.iter().enumerate() {
        let model_number = model_idx + 1;

        let vdw = fitness::vdw_energy(&receptor, ligand);
        let elec = fitness::elec_energy(&receptor, ligand);
        let desolv = fitness::desolv_energy(&receptor, ligand);

        // Calculate clash percentage
        let clash_result = evaluator::calculate_clashes(&receptor, ligand);

        // Calculate AIR energy only if restraints are provided
        let air = match &restraint_pairs {
            Some(pairs) => {
                let restraints = restraints::create_restraints_from_pairs(&receptor, ligand, pairs);
                fitness::air_energy(&restraints, &receptor, ligand)
            }
            None => 0.0,
        };

        let total_score =
            weights.vdw * vdw + weights.elec * elec + weights.desolv * desolv + weights.air * air;

        // Build output line with model number, score, energy terms, then weights
        let mut output = format!(
            "{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            model_number, total_score, vdw, elec, desolv
        );

        if restraint_pairs.is_some() {
            output.push_str(&format!("\t{:.3}", air));
        }

        // Add clash percentage
        output.push_str(&format!("\t{:.2}", clash_result.clash_percentage));

        output.push_str(&format!(
            "\t{}\t{}\t{}\t{}",
            weights.vdw, weights.elec, weights.desolv, weights.air
        ));

        // Calculate DockQ metrics if reference is provided
        if let Some(ref_file) = &reference_file {
            let (_, reference_ligand) = scoring::read_complex(ref_file);
            let evaluator = evaluator::Evaluator::new(receptor.clone(), reference_ligand);
            let metrics = evaluator.calc_metrics(ligand);

            output.push_str(&format!(
                "\t{:.3}\t{:.2}\t{:.2}\t{:.3}",
                metrics.dockq, metrics.rmsd, metrics.irmsd, metrics.fnat
            ));
        }

        println!("{}", output);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants;

    #[test]
    fn test_score_basic() {
        // This test just verifies the function doesn't panic
        // Actual output goes to stdout
        let weights = constants::EnergyWeights::default();

        // We can't easily capture stdout in a unit test without extra crates,
        // so we just verify the function runs without panicking
        score(
            "data/A.pdb".to_string(),
            "data/B.pdb".to_string(),
            None,
            None,
            weights,
        );
    }
}
