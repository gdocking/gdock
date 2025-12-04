pub mod chromosome;
pub mod constants;
pub mod evaluator;
pub mod fitness;
pub mod ga;
pub mod metrics;
pub mod population;
pub mod restraints;
pub mod scoring;
pub mod structure;
pub mod toppar;
pub mod utils;
use rand::rngs::StdRng;
use rand::SeedableRng;

use clap::Command;
use colored::*;
use constants::{
    CONVERGENCE_THRESHOLD, CONVERGENCE_WINDOW, DEFAULT_W_AIR, DEFAULT_W_DESOLV, DEFAULT_W_ELEC,
    DEFAULT_W_VDW, ENABLE_EARLY_STOPPING, MAX_GENERATIONS, POPULATION_SIZE,
};
use indicatif::{ProgressBar, ProgressStyle};
use structure::read_pdb;

fn score(
    receptor_file: String,
    ligand_file: String,
    restraint_pairs: Vec<(i32, i32)>,
    weights: constants::EnergyWeights,
) {
    let receptor = read_pdb(&receptor_file);
    let ligand = read_pdb(&ligand_file);

    let restraints = restraints::create_restraints_from_pairs(&receptor, &ligand, &restraint_pairs);

    let vdw = fitness::vdw_energy(&receptor, &ligand);
    let elec = fitness::elec_energy(&receptor, &ligand);
    let desolv = fitness::desolv_energy(&receptor, &ligand);
    let air = fitness::air_energy(&restraints, &receptor, &ligand);

    let total_score =
        weights.vdw * vdw + weights.elec * elec + weights.desolv * desolv + weights.air * air;

    println!(
        "vdw: {:.3} (w={}) elec: {:.3} (w={}) desolv: {:.3} (w={}) air: {:.3} (w={})",
        vdw, weights.vdw, elec, weights.elec, desolv, weights.desolv, air, weights.air,
    );
    println!("total score: {:.3}", total_score);
}

fn run(
    receptor_file: String,
    ligand_file: String,
    restraint_pairs: Vec<(i32, i32)>,
    reference_file: Option<String>,
    weights: constants::EnergyWeights,
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

    // --------------------------------------------------------------------------------
    // Define number of threads either with the function below or with the
    // environment variable RAYON_NUM_THREADS
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(1)
    //     .build_global()
    //     .unwrap();
    // --------------------------------------------------------------------------------

    let receptor = read_pdb(&receptor_file);
    let ligand = read_pdb(&ligand_file);

    // Create restraints from user-specified residue pairs
    let restraints = restraints::create_restraints_from_pairs(&receptor, &ligand, &restraint_pairs);
    let num_restraints = restraints.len();
    println!(
        "{} Created {} distance restraints\n",
        "✓".green(),
        num_restraints.to_string().cyan()
    );

    // Clone the original molecule for potential RMSD calculations
    let orig = ligand.clone();

    let ligand = utils::position_ligand(&receptor, ligand);

    // Start the evaluator (only if reference is provided)
    let evaluator = if let Some(ref_file) = &reference_file {
        let (_, reference_ligand) = scoring::read_complex(ref_file);
        Some(evaluator::Evaluator::new(
            receptor.clone(),
            reference_ligand,
        ))
    } else {
        None
    };

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

    // Start the GA
    // Create the initial population
    let mut population =
        population::Population::new(Vec::new(), receptor, ligand, orig, restraints, weights);
    let mut rng = StdRng::seed_from_u64(constants::RANDOM_SEED);
    for _ in 0..POPULATION_SIZE {
        let c = chromosome::Chromosome::new(&mut rng);
        population.chromosomes.push(c);
    }

    // Evolve the population
    let mut generation_count = 0;
    let mut best_score_history: Vec<f64> = Vec::new();
    let mut generations_without_improvement = 0;
    let mut last_best_score = f64::MAX;

    while generation_count < MAX_GENERATIONS {
        population.eval_fitness();

        // Calculate metrics for all chromosomes (only if reference is available)
        let metric_vec = evaluator.as_ref().map(|eval| population.eval_metrics(eval));

        // Find the best fitness chromosome and its corresponding RMSD
        let best_fitness_idx = population
            .chromosomes
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.fitness.partial_cmp(&b.fitness).unwrap())
            .map(|(idx, _)| idx)
            .unwrap();

        let best_fitness = population.chromosomes[best_fitness_idx].fitness;
        let best_chr = &population.chromosomes[best_fitness_idx];

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
                "{} Converged at generation {} (no improvement for {} gens)",
                "✓".green(),
                generation_count,
                CONVERGENCE_WINDOW
            ));
            println!();
            break;
        }

        // Calculate population averages
        let mean_fitness = population.get_mean_fitness();
        let mean_rest: f64 = 100.0
            * (population
                .chromosomes
                .iter()
                .map(|c| 1.0 - c.restraint_penalty / num_restraints as f64)
                .sum::<f64>()
                / population.chromosomes.len() as f64);

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

            progress.set_message(format!(
                "DockQ: {:.3} | Score: {:.0}",
                best_metrics.dockq, best_fitness
            ));

            progress.println(format!("  [{}] {} score={:>8.1} dockq={} rest={}% │ {} score={:>8.1} dockq={} rest={}% rmsd={:.2}Å fnat={:.3} irmsd={:.2}Å │ Δ={}%",
                format!("{:>3}", generation_count).bright_black(),
                "📊".bright_blue(),
                mean_fitness,
                format!("{:.3}", mean_dockq).cyan(),
                format!("{:>3.0}", mean_rest).bright_black(),
                "🎯".bright_green(),
                best_fitness,
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
        population = population.evolve(&mut rng);
    }

    // Final evaluation for saving models
    population.eval_fitness();

    let best_fitness_idx = population
        .chromosomes
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.fitness.partial_cmp(&b.fitness).unwrap())
        .map(|(idx, _)| idx)
        .unwrap();

    let final_best_score = &population.chromosomes[best_fitness_idx];

    // Finish progress bar if not already finished
    if !progress.is_finished() {
        progress.finish_and_clear();
    }

    // Save the best models to disk
    println!("\n{}", "💾 Saving Results".bold().cyan());

    // Apply transformations and save best-by-score model
    let best_score_ligand = final_best_score.apply_genes(&ligand_clone);
    let best_score_complex = combine_molecules(&receptor_clone, &best_score_ligand);
    structure::write_pdb(&best_score_complex, &"best_by_score.pdb".to_string());

    // If we have a reference, also save best-by-DockQ model
    if let Some(ref eval) = evaluator {
        let final_metrics = population.eval_metrics(eval);
        let best_dockq_idx = final_metrics
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.dockq.partial_cmp(&b.dockq).unwrap())
            .map(|(idx, _)| idx)
            .unwrap();

        let final_best_dockq = &population.chromosomes[best_dockq_idx];
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

/// Combines receptor and ligand into a single molecule for PDB output
fn combine_molecules(
    receptor: &structure::Molecule,
    ligand: &structure::Molecule,
) -> structure::Molecule {
    let mut combined = structure::Molecule::new();

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

fn main() {
    const VERSION: &str = env!("CARGO_PKG_VERSION");
    let matches = Command::new("gdock")
        .version(VERSION)
        .author("pending")
        .about("Fast information-driven protein-protein docking using genetic algorithms")
        .arg(
            clap::Arg::new("receptor")
                .long("receptor")
                .short('r')
                .value_name("FILE")
                .help("Receptor PDB file")
                .required(true),
        )
        .arg(
            clap::Arg::new("ligand")
                .long("ligand")
                .short('l')
                .value_name("FILE")
                .help("Ligand PDB file")
                .required(true),
        )
        .arg(
            clap::Arg::new("restraints")
                .long("restraints")
                .value_name("PAIRS")
                .help("Comma-separated restraint pairs receptor:ligand (e.g., 10:45,15:50,20:55)")
                .required(true),
        )
        .arg(
            clap::Arg::new("reference")
                .long("reference")
                .value_name("FILE")
                .help("Reference PDB file for DockQ calculation (optional)"),
        )
        .arg(
            clap::Arg::new("score")
                .long("score")
                .help("Score mode: only calculate energy without running GA")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            clap::Arg::new("w_vdw")
                .long("w_vdw")
                .value_name("WEIGHT")
                .help("Weight for VDW energy term (default: 1.0)")
                .value_parser(clap::value_parser!(f64)),
        )
        .arg(
            clap::Arg::new("w_elec")
                .long("w_elec")
                .value_name("WEIGHT")
                .help("Weight for electrostatic energy term (default: 0.5)")
                .value_parser(clap::value_parser!(f64)),
        )
        .arg(
            clap::Arg::new("w_desolv")
                .long("w_desolv")
                .value_name("WEIGHT")
                .help("Weight for desolvation energy term (default: 0.5)")
                .value_parser(clap::value_parser!(f64)),
        )
        .arg(
            clap::Arg::new("w_air")
                .long("w_air")
                .value_name("WEIGHT")
                .help("Weight for AIR restraint energy term (default: 100.0)")
                .value_parser(clap::value_parser!(f64)),
        )
        .get_matches();

    // Parse common arguments
    let receptor_file = matches.get_one::<String>("receptor").unwrap().clone();
    let ligand_file = matches.get_one::<String>("ligand").unwrap().clone();
    let reference_file = matches.get_one::<String>("reference").cloned();
    let score_only = matches.get_flag("score");

    // Parse restraint pairs (format: "rec1:lig1,rec2:lig2,...")
    let restraint_pairs: Vec<(i32, i32)> = matches
        .get_one::<String>("restraints")
        .unwrap()
        .split(',')
        .map(|pair| {
            let parts: Vec<&str> = pair.trim().split(':').collect();
            if parts.len() != 2 {
                panic!(
                    "Invalid restraint format: '{}'. Expected format: 'receptor:ligand'",
                    pair
                );
            }
            let rec = parts[0]
                .trim()
                .parse::<i32>()
                .unwrap_or_else(|_| panic!("Invalid receptor residue number: '{}'", parts[0]));
            let lig = parts[1]
                .trim()
                .parse::<i32>()
                .unwrap_or_else(|_| panic!("Invalid ligand residue number: '{}'", parts[1]));
            (rec, lig)
        })
        .collect();

    // Parse energy weights
    let weights = constants::EnergyWeights::new(
        matches
            .get_one::<f64>("w_vdw")
            .copied()
            .unwrap_or(DEFAULT_W_VDW),
        matches
            .get_one::<f64>("w_elec")
            .copied()
            .unwrap_or(DEFAULT_W_ELEC),
        matches
            .get_one::<f64>("w_desolv")
            .copied()
            .unwrap_or(DEFAULT_W_DESOLV),
        matches
            .get_one::<f64>("w_air")
            .copied()
            .unwrap_or(DEFAULT_W_AIR),
    );

    // Execute based on mode
    if score_only {
        score(receptor_file, ligand_file, restraint_pairs, weights);
    } else {
        run(
            receptor_file,
            ligand_file,
            restraint_pairs,
            reference_file,
            weights,
        );
    }
}
