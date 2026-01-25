use colored::*;
use indicatif::{ProgressBar, ProgressStyle};
use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::chromosome;
use crate::constants::{
    self, EnergyWeights, CONVERGENCE_THRESHOLD, CONVERGENCE_WINDOW, ENABLE_EARLY_STOPPING,
    MAX_GENERATIONS, POPULATION_SIZE,
};
use crate::evaluator;
use crate::population;
use crate::restraints;
use crate::scoring;
use crate::structure::{self, read_pdb, Molecule};
use crate::utils;

/// Combines receptor and ligand into a single molecule for PDB output
pub fn combine_molecules(receptor: &Molecule, ligand: &Molecule) -> Molecule {
    let mut combined = Molecule::new();

    // Add receptor atoms
    for atom in &receptor.0 {
        combined.0.push(atom.clone());
    }

    // Add ligand atoms
    for atom in &ligand.0 {
        combined.0.push(atom.clone());
    }

    combined
}

/// Run the genetic algorithm docking
pub fn run(
    receptor_file: String,
    ligand_file: String,
    restraint_pairs: Vec<(i32, i32)>,
    reference_file: Option<String>,
    weights: EnergyWeights,
    debug_mode: bool,
) {
    const VERSION: &str = env!("CARGO_PKG_VERSION");
    println!(
        "\n{} {}",
        "🧬 GDock".bold().cyan(),
        format!("v{}", VERSION).bright_black()
    );
    println!(
        "{}",
        "   Protein-Protein Docking with Genetic Algorithm".bright_black()
    );
    if debug_mode {
        println!(
            "{}",
            "   ⚠️  DEBUG MODE: Using DockQ as fitness function"
                .yellow()
                .bold()
        );
    }
    println!(
        "{}",
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━".bright_black()
    );
    println!("{}", "📁 Input Files".bold());
    println!("  {}  {}", "Receptor: ".green(), receptor_file);
    println!("  {}    {}", "Ligand: ".green(), ligand_file);
    if let Some(ref_file) = &reference_file {
        println!(
            "  {} {} {}",
            "Reference:".green(),
            ref_file,
            "(DockQ mode)".bright_black()
        );
    } else {
        println!(
            "  {} {}",
            "Reference:".green(),
            "None (score-only mode)".yellow()
        );
    }
    println!(
        "\n{} {} pairs",
        "🎯 Restraints:".bold(),
        restraint_pairs.len().to_string().cyan()
    );
    for (rec, lig) in restraint_pairs.iter() {
        println!("  {} {}:{}", "•".bright_blue(), rec, lig);
    }
    println!("\n{}", "⚙️  Energy Weights".bold());
    println!(
        "  {}={:.2} │ {}={:.2} │ {}={:.2} │ {}={:.2}",
        "VDW".green(),
        weights.vdw,
        "Elec".green(),
        weights.elec,
        "Desolv".green(),
        weights.desolv,
        "AIR".green(),
        weights.air
    );
    println!(
        "{}\n",
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━".bright_black()
    );

    let receptor_model = read_pdb(&receptor_file);
    let ligand_model = read_pdb(&ligand_file);

    // For docking mode, use the first model
    let receptor = receptor_model.0[0].clone();
    let ligand = ligand_model.0[0].clone();

    // Create restraints from user-specified residue pairs
    let restraints_list =
        restraints::create_restraints_from_pairs(&receptor, &ligand, &restraint_pairs);
    let num_restraints = restraints_list.len();
    println!(
        "{} Created {} distance restraints\n",
        "✓".green(),
        num_restraints.to_string().cyan()
    );

    // Clone the original molecule for potential RMSD calculations
    let orig = ligand.clone();

    let ligand = utils::position_ligand(&receptor, ligand);

    // Start the evaluator (only if reference is provided)
    let eval = if let Some(ref_file) = &reference_file {
        let (_, reference_ligand) = scoring::read_complex(ref_file);
        Some(evaluator::Evaluator::new(
            receptor.clone(),
            reference_ligand,
        ))
    } else {
        None
    };

    // In debug mode, we use the evaluator as the fitness function
    let debug_evaluator = if debug_mode { eval.clone() } else { None };

    // Clone molecules before moving them into population (needed for PDB output later)
    let receptor_clone = receptor.clone();
    let ligand_clone = ligand.clone();

    // Create progress bar
    let progress = ProgressBar::new(MAX_GENERATIONS);
    progress.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{bar:40.cyan/blue}] {pos}/{len} gens | {msg}")
            .unwrap()
            .progress_chars("█▓░"),
    );

    // Create the initial population
    let mut pop = population::Population::new(
        Vec::new(),
        receptor,
        ligand,
        orig,
        restraints_list,
        weights,
        debug_evaluator,
    );
    let mut rng = StdRng::seed_from_u64(constants::RANDOM_SEED);
    for _ in 0..POPULATION_SIZE {
        let c = chromosome::Chromosome::new(&mut rng);
        pop.chromosomes.push(c);
    }

    // Evolve the population
    let mut generation_count = 0;
    let mut best_score_history: Vec<f64> = Vec::new();
    let mut generations_without_improvement = 0;
    let mut last_best_score = f64::MAX;

    println!(
        "{} Starting evolution for {} generations",
        "🧬".bold(),
        MAX_GENERATIONS
    );
    while generation_count < MAX_GENERATIONS {
        pop.eval_fitness();

        // Calculate metrics for all chromosomes (only if reference is available)
        let metric_vec = eval.as_ref().map(|e| pop.eval_metrics(e));

        // Find the best fitness chromosome
        let best_fitness_idx = pop
            .chromosomes
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.fitness.partial_cmp(&b.fitness).unwrap())
            .map(|(idx, _)| idx)
            .unwrap();

        let best_fitness = pop.chromosomes[best_fitness_idx].fitness;
        let best_chr = &pop.chromosomes[best_fitness_idx];

        // Track convergence
        best_score_history.push(best_fitness);
        if generation_count > 0 {
            let improvement = (last_best_score - best_fitness) / last_best_score.abs();
            if improvement < CONVERGENCE_THRESHOLD {
                generations_without_improvement += 1;
            } else {
                generations_without_improvement = 0;
            }
        }
        last_best_score = best_fitness;

        // Check for early stopping
        if ENABLE_EARLY_STOPPING && generations_without_improvement >= CONVERGENCE_WINDOW {
            progress.finish_with_message(format!(
                "{} Converged at generation {} (no improvement larger than {}% for {} gens)",
                "✓".green(),
                generation_count,
                CONVERGENCE_THRESHOLD * 100.0,
                CONVERGENCE_WINDOW
            ));
            println!();
            break;
        }

        // Calculate population averages
        let mean_fitness = pop.get_mean_fitness();
        let mean_rest: f64 = 100.0
            * (pop
                .chromosomes
                .iter()
                .map(|c| 1.0 - c.restraint_penalty / num_restraints as f64)
                .sum::<f64>()
                / pop.chromosomes.len() as f64);

        let best_rest = 100.0 * (1.0 - best_chr.restraint_penalty / num_restraints as f64);

        // Calculate improvement metrics
        let improvement_since_last = if generation_count > 0 {
            let prev = best_score_history[generation_count as usize - 1];
            ((prev - best_fitness) / prev.abs()) * 100.0
        } else {
            0.0
        };

        // Update progress bar and print output
        progress.set_position(generation_count);

        if let Some(ref metrics) = metric_vec {
            let mean_dockq: f64 =
                metrics.iter().map(|m| m.dockq).sum::<f64>() / metrics.len() as f64;
            let best_metrics = &metrics[best_fitness_idx];

            let dockq_color = if best_metrics.dockq >= 0.8 {
                "green"
            } else if best_metrics.dockq >= 0.5 {
                "yellow"
            } else if best_metrics.dockq >= 0.23 {
                "bright_yellow"
            } else {
                "red"
            };

            let score_label = if debug_mode { "DockQ" } else { "Score" };
            let score_value = if debug_mode {
                -best_fitness
            } else {
                best_fitness
            };
            progress.set_message(format!(
                "DockQ: {:.3} | {}: {:.3}",
                best_metrics.dockq, score_label, score_value
            ));

            let mean_score_display = if debug_mode {
                -mean_fitness
            } else {
                mean_fitness
            };
            let best_score_display = if debug_mode {
                -best_fitness
            } else {
                best_fitness
            };
            progress.println(format!("  [{}] {} score={:>8.3} dockq={} rest={}% │ {} score={:>8.3} dockq={} rest={}% rmsd={:.2}Å fnat={:.3} irmsd={:.2}Å │ Δ={}%",
                format!("{:>3}", generation_count).bright_black(),
                "📊".bright_blue(),
                mean_score_display,
                format!("{:.3}", mean_dockq).cyan(),
                format!("{:>3.0}", mean_rest).bright_black(),
                "🎯".bright_green(),
                best_score_display,
                match dockq_color {
                    "green" => format!("{:.3}", best_metrics.dockq).green(),
                    "yellow" => format!("{:.3}", best_metrics.dockq).yellow(),
                    "bright_yellow" => format!("{:.3}", best_metrics.dockq).bright_yellow(),
                    _ => format!("{:.3}", best_metrics.dockq).red()
                },
                format!("{:>3.0}", best_rest).bright_black(),
                best_metrics.rmsd,
                best_metrics.fnat,
                best_metrics.irmsd,
                if improvement_since_last > 0.0 {
                    format!("{:>+5.2}", improvement_since_last).green()
                } else {
                    format!("{:>+5.2}", improvement_since_last).bright_black()
                }
            ));
        } else {
            progress.set_message(format!("Score: {:.0}", best_fitness));
            progress.println(format!(
                "  [{}] {} score={:>8.1} rest={}% │ {} score={:>8.1} rest={}% │ Δ={}%",
                format!("{:>3}", generation_count).bright_black(),
                "📊".bright_blue(),
                mean_fitness,
                format!("{:>3.0}", mean_rest).bright_black(),
                "🎯".bright_green(),
                best_fitness,
                format!("{:>3.0}", best_rest).bright_black(),
                if improvement_since_last > 0.0 {
                    format!("{:>+5.2}", improvement_since_last).green()
                } else {
                    format!("{:>+5.2}", improvement_since_last).bright_black()
                }
            ));
        }

        generation_count += 1;
        pop = pop.evolve(&mut rng);
    }

    // Final evaluation for saving models
    pop.eval_fitness();

    let best_fitness_idx = pop
        .chromosomes
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.fitness.partial_cmp(&b.fitness).unwrap())
        .map(|(idx, _)| idx)
        .unwrap();

    let final_best_score = &pop.chromosomes[best_fitness_idx];

    // Finish progress bar if not already finished
    if !progress.is_finished() {
        progress.finish();
    }

    // Save the best models to disk
    println!("\n{}", "💾 Saving Results".bold().cyan());

    // Display final metrics
    if let Some(ref e) = eval {
        let final_metrics = pop.eval_metrics(e);
        let best_dockq_idx = final_metrics
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.dockq.partial_cmp(&b.dockq).unwrap())
            .map(|(idx, _)| idx)
            .unwrap();
        let best_score_metrics = &final_metrics[best_fitness_idx];
        let best_dockq_metrics = &final_metrics[best_dockq_idx];

        println!("\n{}", "📊 Final Metrics".bold().cyan());
        println!(
            "  {} DockQ={:.3} RMSD={:.2}Å iRMSD={:.2}Å FNAT={:.3}",
            "Best by score:".green(),
            best_score_metrics.dockq,
            best_score_metrics.rmsd,
            best_score_metrics.irmsd,
            best_score_metrics.fnat
        );
        println!(
            "  {} DockQ={:.3} RMSD={:.2}Å iRMSD={:.2}Å FNAT={:.3}",
            "Best by DockQ:".green(),
            best_dockq_metrics.dockq,
            best_dockq_metrics.rmsd,
            best_dockq_metrics.irmsd,
            best_dockq_metrics.fnat
        );
    }

    // Apply transformations and save best-by-score model
    let best_score_ligand = final_best_score.apply_genes(&ligand_clone);
    let best_score_complex = combine_molecules(&receptor_clone, &best_score_ligand);
    structure::write_pdb(&best_score_complex, &"best_by_score.pdb".to_string());

    // If we have a reference, also save best-by-DockQ model
    if let Some(ref e) = eval {
        let final_metrics = pop.eval_metrics(e);
        let best_dockq_idx = final_metrics
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.dockq.partial_cmp(&b.dockq).unwrap())
            .map(|(idx, _)| idx)
            .unwrap();

        let final_best_dockq = &pop.chromosomes[best_dockq_idx];
        let best_dockq_ligand = final_best_dockq.apply_genes(&ligand_clone);
        let best_dockq_complex = combine_molecules(&receptor_clone, &best_dockq_ligand);
        structure::write_pdb(&best_dockq_complex, &"best_by_dockq.pdb".to_string());

        println!("  {} {}", "✓".green(), "best_by_score.pdb".bright_white());
        println!("  {} {}", "✓".green(), "best_by_dockq.pdb".bright_white());
    } else {
        println!("  {} {}", "✓".green(), "best_by_score.pdb".bright_white());
    }
    println!("\n{}", "✨ Done!".bold().green());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_combine_molecules() {
        let receptor_model = read_pdb(&"data/A.pdb".to_string());
        let ligand_model = read_pdb(&"data/B.pdb".to_string());

        let receptor = &receptor_model.0[0];
        let ligand = &ligand_model.0[0];

        let combined = combine_molecules(receptor, ligand);

        assert_eq!(
            combined.0.len(),
            receptor.0.len() + ligand.0.len(),
            "Combined molecule should have all atoms from both"
        );
    }
}
