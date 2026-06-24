use colored::*;
use indicatif::{ProgressBar, ProgressStyle};
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::cmp::Ordering;
use std::fs;
use std::io::Write;
use std::path::PathBuf;

use crate::chromosome;
use crate::constants::{
    self, EnergyWeights, CONVERGENCE_THRESHOLD, CONVERGENCE_WINDOW, HALL_OF_FAME_MAX_SIZE,
    MAX_GENERATIONS, POPULATION_SIZE,
};
use crate::evaluator;
use crate::hall_of_fame::{HallOfFame, HallOfFameEntry};
use crate::population;
use crate::restraints;
use crate::runner::{run_ga, select_models};
use crate::scoring;
use crate::structure::{self, read_pdb, Molecule};
use crate::utils;

/// Configuration for a docking run
pub struct RunConfig {
    pub receptor_file: String,
    pub ligand_file: String,
    pub restraint_pairs: Vec<(i32, i32)>,
    pub reference_file: Option<String>,
    pub weights: EnergyWeights,
    pub debug_mode: bool,
    pub output_dir: Option<String>,
    pub no_clustering: bool,
    pub sampling: Option<usize>,
}

/// Re-exported for use by tests and other modules that imported from here.
pub use crate::structure::combine_molecules;

// ============================================================================

/// Run the genetic algorithm docking
pub fn run(config: RunConfig) {
    let RunConfig {
        receptor_file,
        ligand_file,
        restraint_pairs,
        reference_file,
        weights,
        debug_mode,
        output_dir,
        no_clustering,
        sampling,
    } = config;
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

    // Create bidirectional ambiguous restraints (HADDOCK-style: anchor on each chain)
    let restraints_list =
        restraints::create_ambiguous_restraints_from_pairs(&receptor, &ligand, &restraint_pairs);
    let num_restraints = restraints_list.len() / 2;
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
    let mut best_score_history: Vec<f64> = Vec::new();

    println!(
        "{} Starting evolution for {} generations",
        "🧬".bold(),
        MAX_GENERATIONS
    );

    let hof_capacity = sampling.unwrap_or(HALL_OF_FAME_MAX_SIZE);
    let ga_result = run_ga(pop, &mut rng, MAX_GENERATIONS, hof_capacity, |gen, pop| {
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

        best_score_history.push(best_fitness);

        let improvement_since_last = if gen > 0 {
            let prev = best_score_history[gen as usize - 1];
            if prev.abs() < f64::EPSILON {
                0.0
            } else {
                ((prev - best_fitness) / prev.abs()) * 100.0
            }
        } else {
            0.0
        };

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

        // Update progress bar and print output
        progress.set_position(gen);

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
                format!("{:>3}", gen).bright_black(),
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
                format!("{:>3}", gen).bright_black(),
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
    });

    // Finish progress bar
    if ga_result.converged_early {
        progress.finish_with_message(format!(
            "{} Converged at generation {} (no improvement larger than {}% for {} gens)",
            "✓".green(),
            ga_result.generations_run,
            CONVERGENCE_THRESHOLD * 100.0,
            CONVERGENCE_WINDOW
        ));
        println!();
    } else {
        progress.finish();
    }

    let hall_of_fame = ga_result.hall_of_fame;
    let pop = ga_result.final_population;

    // Determine output directory
    let out_dir = match &output_dir {
        Some(dir) => {
            let path = PathBuf::from(dir);
            fs::create_dir_all(&path).expect("Failed to create output directory");
            path
        }
        None => PathBuf::from("."),
    };

    if no_clustering {
        // =====================================================================
        // No clustering: Output best_by_score and best_by_dockq (old behavior)
        // =====================================================================

        let best_fitness_idx = pop
            .chromosomes
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.fitness.partial_cmp(&b.fitness).unwrap())
            .map(|(idx, _)| idx)
            .unwrap();

        let final_best_score = &pop.chromosomes[best_fitness_idx];

        println!("\n{}", "💾 Saving Results".bold().cyan());

        // Save best-by-score model
        let best_score_ligand = final_best_score.apply_genes(&ligand_clone);
        let best_score_complex = combine_molecules(&receptor_clone, &best_score_ligand);
        let best_score_path = out_dir.join("best_by_score.pdb");
        structure::write_pdb(
            &best_score_complex,
            best_score_path.to_string_lossy().as_ref(),
        );

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
            let best_dockq_path = out_dir.join("best_by_dockq.pdb");
            structure::write_pdb(
                &best_dockq_complex,
                best_dockq_path.to_string_lossy().as_ref(),
            );

            // Write metrics.tsv
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

            let metrics_path = out_dir.join("metrics.tsv");
            let mut metrics_file =
                fs::File::create(&metrics_path).expect("Failed to create metrics file");
            writeln!(metrics_file, "model\tdockq\trmsd\tirmsd\tfnat\tscore").unwrap();
            writeln!(
                metrics_file,
                "best_by_score\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                best_score_metrics.dockq,
                best_score_metrics.rmsd,
                best_score_metrics.irmsd,
                best_score_metrics.fnat,
                final_best_score.fitness
            )
            .unwrap();
            writeln!(
                metrics_file,
                "best_by_dockq\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                best_dockq_metrics.dockq,
                best_dockq_metrics.rmsd,
                best_dockq_metrics.irmsd,
                best_dockq_metrics.fnat,
                final_best_dockq.fitness
            )
            .unwrap();

            println!("  {} {}", "✓".green(), best_score_path.display());
            println!("  {} {}", "✓".green(), best_dockq_path.display());
            println!("  {} {}", "✓".green(), metrics_path.display());
        } else {
            println!("  {} {}", "✓".green(), best_score_path.display());
        }
    } else {
        // =====================================================================
        // Clustering: Select diverse representative structures
        // =====================================================================

        // Report Hall of Fame status
        println!(
            "\n{} Collected {} diverse structures in Hall of Fame",
            "📦".bold(),
            hall_of_fame.len().to_string().cyan()
        );

        println!(
            "\n{}",
            "🔬 Clustering Hall of Fame structures".bold().cyan()
        );

        let hof_entries = hall_of_fame.entries();
        let selected = select_models(hof_entries, &receptor_clone, &ligand_clone);

        // Save output models and metrics
        println!("\n{}", "💾 Saving Results".bold().cyan());

        let metrics_path = out_dir.join("metrics.tsv");
        let mut metrics_file =
            fs::File::create(&metrics_path).expect("Failed to create metrics file");

        if eval.is_some() {
            writeln!(
                metrics_file,
                "model\tcluster_size\tscore\tdockq\trmsd\tirmsd\tfnat"
            )
            .unwrap();
        } else {
            writeln!(metrics_file, "model\tcluster_size\tscore").unwrap();
        }

        println!("\n{}", "📊 Output Models (FCC Clustered)".bold().cyan());

        for (model_num, (hof_idx, cluster_size)) in selected.clustered.iter().enumerate() {
            let entry = &hof_entries[*hof_idx];
            let model_name = format!("model_{}", model_num + 1);

            let ligand = ligand_clone
                .clone()
                .rotate(entry.genes[0], entry.genes[1], entry.genes[2])
                .displace(entry.genes[3], entry.genes[4], entry.genes[5]);
            let complex = combine_molecules(&receptor_clone, &ligand);

            let pdb_path = out_dir.join(format!("{}.pdb", model_name));
            structure::write_pdb(&complex, pdb_path.to_string_lossy().as_ref());

            if let Some(ref e) = eval {
                let metrics = e.calc_metrics(&ligand);

                writeln!(
                    metrics_file,
                    "{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                    model_name,
                    cluster_size,
                    entry.fitness,
                    metrics.dockq,
                    metrics.rmsd,
                    metrics.irmsd,
                    metrics.fnat
                )
                .unwrap();

                let dockq_str = if metrics.dockq >= 0.80 {
                    format!("{:.3}", metrics.dockq).green()
                } else if metrics.dockq >= 0.49 {
                    format!("{:.3}", metrics.dockq).yellow()
                } else if metrics.dockq >= 0.23 {
                    format!("{:.3}", metrics.dockq).bright_yellow()
                } else {
                    format!("{:.3}", metrics.dockq).red()
                };

                println!(
                    "  {}: cluster={} score={:.1} DockQ={}",
                    model_name.green(),
                    format!("{:>3}", cluster_size).cyan(),
                    entry.fitness,
                    dockq_str
                );
            } else {
                writeln!(
                    metrics_file,
                    "{}\t{}\t{:.4}",
                    model_name, cluster_size, entry.fitness
                )
                .unwrap();

                println!(
                    "  {}: cluster={} score={:.1}",
                    model_name.green(),
                    format!("{:>3}", cluster_size).cyan(),
                    entry.fitness
                );
            }

            println!("    {} {}", "✓".bright_black(), pdb_path.display());
        }

        // =====================================================================
        // Also output top 5 by score (ranked_*.pdb)
        // =====================================================================

        println!("\n{}", "📊 Output Models (Ranked by Score)".bold().cyan());

        for (rank, hof_idx) in selected.ranked.iter().enumerate() {
            let entry = &hof_entries[*hof_idx];
            let model_name = format!("ranked_{}", rank + 1);

            let ligand = ligand_clone
                .clone()
                .rotate(entry.genes[0], entry.genes[1], entry.genes[2])
                .displace(entry.genes[3], entry.genes[4], entry.genes[5]);
            let complex = combine_molecules(&receptor_clone, &ligand);

            let pdb_path = out_dir.join(format!("{}.pdb", model_name));
            structure::write_pdb(&complex, pdb_path.to_string_lossy().as_ref());

            if let Some(ref e) = eval {
                let metrics = e.calc_metrics(&ligand);

                writeln!(
                    metrics_file,
                    "{}\t-\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                    model_name,
                    entry.fitness,
                    metrics.dockq,
                    metrics.rmsd,
                    metrics.irmsd,
                    metrics.fnat
                )
                .unwrap();

                let dockq_str = if metrics.dockq >= 0.80 {
                    format!("{:.3}", metrics.dockq).green()
                } else if metrics.dockq >= 0.49 {
                    format!("{:.3}", metrics.dockq).yellow()
                } else if metrics.dockq >= 0.23 {
                    format!("{:.3}", metrics.dockq).bright_yellow()
                } else {
                    format!("{:.3}", metrics.dockq).red()
                };

                println!(
                    "  {}: score={:.1} DockQ={}",
                    model_name.green(),
                    entry.fitness,
                    dockq_str
                );
            } else {
                writeln!(metrics_file, "{}\t-\t{:.4}", model_name, entry.fitness).unwrap();

                println!("  {}: score={:.1}", model_name.green(), entry.fitness);
            }

            println!("    {} {}", "✓".bright_black(), pdb_path.display());
        }

        println!("    {} {}", "✓".bright_black(), metrics_path.display());
    }

    if sampling.is_some() {
        write_sampling_output(&hall_of_fame, &out_dir, &receptor_clone, &ligand_clone);
    }

    println!("\n{}", "✨ Done!".bold().green());
}

fn write_sampling_output(
    hall_of_fame: &HallOfFame,
    out_dir: &std::path::Path,
    receptor: &Molecule,
    ligand: &Molecule,
) {
    let sampling_dir = out_dir.join("sampling");
    fs::create_dir_all(&sampling_dir).expect("Failed to create sampling directory");

    let mut sorted: Vec<&HallOfFameEntry> = hall_of_fame.entries().iter().collect();
    sorted.sort_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap_or(Ordering::Equal));

    let tsv_path = sampling_dir.join("sampling.tsv");
    let mut tsv = fs::File::create(&tsv_path).expect("Failed to create sampling.tsv");
    writeln!(tsv, "model\tscore\tvdw\telec\tdesolv\tair").unwrap();

    println!("\n{}", "📦 Sampling output".bold().cyan());

    for (rank, entry) in sorted.iter().enumerate() {
        let model_name = format!("gdock_{}", rank + 1);
        let docked_ligand = ligand
            .clone()
            .rotate(entry.genes[0], entry.genes[1], entry.genes[2])
            .displace(entry.genes[3], entry.genes[4], entry.genes[5]);
        let complex = combine_molecules(receptor, &docked_ligand);
        let pdb_path = sampling_dir.join(format!("{}.pdb", model_name));
        structure::write_pdb(&complex, pdb_path.to_string_lossy().as_ref());
        writeln!(
            tsv,
            "{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            model_name, entry.fitness, entry.vdw, entry.elec, entry.desolv, entry.air
        )
        .unwrap();
    }

    println!(
        "  {} {} structures written to {}",
        "✓".green(),
        sorted.len(),
        sampling_dir.display()
    );
    println!("  {} {}", "✓".green(), tsv_path.display());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hall_of_fame::HallOfFame;
    use std::f64::consts::PI;

    #[test]
    fn test_combine_molecules() {
        let receptor_model = read_pdb("data/2oob_A.pdb");
        let ligand_model = read_pdb("data/2oob_B.pdb");

        let receptor = &receptor_model.0[0];
        let ligand = &ligand_model.0[0];

        let combined = combine_molecules(receptor, ligand);

        assert_eq!(
            combined.0.len(),
            receptor.0.len() + ligand.0.len(),
            "Combined molecule should have all atoms from both"
        );
    }

    #[test]
    fn test_write_sampling_output_creates_files() {
        let receptor_model = read_pdb("data/2oob_A.pdb");
        let ligand_model = read_pdb("data/2oob_B.pdb");
        let receptor = receptor_model.0[0].clone();
        let ligand = ligand_model.0[0].clone();

        let mut hof = HallOfFame::new();
        hof.try_add(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], -100.0, 1.0, 2.0, 3.0, 4.0);
        hof.try_add(&[PI, PI, PI, 10.0, 10.0, 10.0], -50.0, 1.0, 2.0, 3.0, 4.0);

        let tmp = crate::utils::get_unique_tempdir();
        write_sampling_output(&hof, tmp.as_path(), &receptor, &ligand);

        let sampling_dir = tmp.as_path().join("sampling");
        assert!(
            sampling_dir.exists(),
            "sampling/ directory should be created"
        );
        assert!(
            sampling_dir.join("gdock_1.pdb").exists(),
            "gdock_1.pdb should exist"
        );
        assert!(
            sampling_dir.join("gdock_2.pdb").exists(),
            "gdock_2.pdb should exist"
        );
        assert!(
            sampling_dir.join("sampling.tsv").exists(),
            "sampling.tsv should exist"
        );
    }

    #[test]
    fn test_write_sampling_output_sorted_by_fitness() {
        let receptor_model = read_pdb(&"data/2oob_A.pdb".to_string());
        let ligand_model = read_pdb(&"data/2oob_B.pdb".to_string());
        let receptor = receptor_model.0[0].clone();
        let ligand = ligand_model.0[0].clone();

        // Add entries out of fitness order — worst first
        let mut hof = HallOfFame::new();
        hof.try_add(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], -50.0, 0.0, 0.0, 0.0, 0.0);
        hof.try_add(&[PI, PI, PI, 10.0, 10.0, 10.0], -100.0, 0.0, 0.0, 0.0, 0.0);

        let tmp = crate::utils::get_unique_tempdir();
        write_sampling_output(&hof, tmp.as_path(), &receptor, &ligand);

        let tsv = std::fs::read_to_string(tmp.as_path().join("sampling/sampling.tsv")).unwrap();
        let mut lines = tsv.lines();
        assert_eq!(
            lines.next().unwrap(),
            "model\tscore\tvdw\telec\tdesolv\tair"
        );

        // gdock_1 should be the best (lowest) score
        let first_row = lines.next().unwrap();
        assert!(
            first_row.starts_with("gdock_1"),
            "first row should be gdock_1"
        );
        let score: f64 = first_row.split('\t').nth(1).unwrap().parse().unwrap();
        assert!(score < -50.0, "gdock_1 should have the best (lowest) score");
    }

    #[test]
    fn test_write_sampling_output_tsv_columns() {
        let receptor_model = read_pdb("data/2oob_A.pdb");
        let ligand_model = read_pdb("data/2oob_B.pdb");
        let receptor = receptor_model.0[0].clone();
        let ligand = ligand_model.0[0].clone();

        let mut hof = HallOfFame::new();
        hof.try_add(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], -99.0, 1.1, 2.2, 3.3, 4.4);

        let tmp = crate::utils::get_unique_tempdir();
        write_sampling_output(&hof, tmp.as_path(), &receptor, &ligand);

        let tsv = std::fs::read_to_string(tmp.as_path().join("sampling/sampling.tsv")).unwrap();
        let mut lines = tsv.lines();
        assert_eq!(
            lines.next().unwrap(),
            "model\tscore\tvdw\telec\tdesolv\tair"
        );

        let data_row = lines.next().unwrap();
        let cols: Vec<&str> = data_row.split('\t').collect();
        assert_eq!(cols.len(), 6, "data rows should have 6 columns");
        assert_eq!(cols[0], "gdock_1");
        assert!((cols[1].parse::<f64>().unwrap() - (-99.0)).abs() < 0.001);
        assert!((cols[2].parse::<f64>().unwrap() - 1.1).abs() < 0.001);
        assert!((cols[3].parse::<f64>().unwrap() - 2.2).abs() < 0.001);
        assert!((cols[4].parse::<f64>().unwrap() - 3.3).abs() < 0.001);
        assert!((cols[5].parse::<f64>().unwrap() - 4.4).abs() < 0.001);
    }
}
