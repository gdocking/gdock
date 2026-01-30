use colored::*;
use indicatif::{ProgressBar, ProgressStyle};
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::fs;
use std::io::Write;
use std::path::PathBuf;

use crate::chromosome;
use crate::clustering::{self, ClusteringConfig};
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
}

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

// ============================================================================
// Hall of Fame - tracks diverse solutions throughout evolution
// ============================================================================

/// Default parameters for Hall of Fame
const HALL_OF_FAME_MAX_SIZE: usize = 500;
const HALL_OF_FAME_TOP_K: usize = 10;
const UNIQUENESS_ROTATION_THRESHOLD: f64 = 0.2; // ~11 degrees
const UNIQUENESS_TRANSLATION_THRESHOLD: f64 = 2.0; // 2 Å

/// Number of output models after clustering
const NUM_OUTPUT_MODELS: usize = 5;

/// Entry in the Hall of Fame storing a solution's genes and metadata
#[derive(Debug, Clone)]
pub struct HallOfFameEntry {
    /// Chromosome genes: [alpha, beta, gamma, dx, dy, dz]
    pub genes: [f64; 6],
    /// Fitness score at time of capture
    pub fitness: f64,
    /// Generation when captured
    pub generation: u32,
}

/// Tracks diverse high-quality solutions throughout evolution
#[derive(Debug)]
pub struct HallOfFame {
    entries: Vec<HallOfFameEntry>,
    max_size: usize,
}

impl HallOfFame {
    /// Create a new Hall of Fame with default max size
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            max_size: HALL_OF_FAME_MAX_SIZE,
        }
    }

    /// Create a new Hall of Fame with custom max size
    pub fn with_capacity(max_size: usize) -> Self {
        Self {
            entries: Vec::new(),
            max_size,
        }
    }

    /// Try to add a new entry if it's unique enough
    /// Returns true if the entry was added
    pub fn try_add(&mut self, genes: &[f64], fitness: f64, generation: u32) -> bool {
        if genes.len() != 6 {
            return false;
        }

        let new_genes: [f64; 6] = [genes[0], genes[1], genes[2], genes[3], genes[4], genes[5]];

        // Check uniqueness against existing entries
        if !self.is_unique(&new_genes) {
            return false;
        }

        // Add the entry
        self.entries.push(HallOfFameEntry {
            genes: new_genes,
            fitness,
            generation,
        });

        // If over capacity, remove worst entries
        if self.entries.len() > self.max_size {
            self.prune();
        }

        true
    }

    /// Check if genes are unique compared to existing entries
    fn is_unique(&self, new_genes: &[f64; 6]) -> bool {
        for entry in &self.entries {
            if Self::genes_are_similar(new_genes, &entry.genes) {
                return false;
            }
        }
        true
    }

    /// Check if two gene sets are similar (potential duplicates)
    fn genes_are_similar(a: &[f64; 6], b: &[f64; 6]) -> bool {
        // Check rotations (with wrap-around handling)
        for i in 0..3 {
            let diff = Self::angular_difference(a[i], b[i]);
            if diff > UNIQUENESS_ROTATION_THRESHOLD {
                return false;
            }
        }

        // Check translations
        for i in 3..6 {
            let diff = (a[i] - b[i]).abs();
            if diff > UNIQUENESS_TRANSLATION_THRESHOLD {
                return false;
            }
        }

        true
    }

    /// Calculate minimum angular difference (handles wrap-around at 2π)
    fn angular_difference(a: f64, b: f64) -> f64 {
        use std::f64::consts::PI;
        let diff = (a - b).abs();
        if diff > PI {
            2.0 * PI - diff
        } else {
            diff
        }
    }

    /// Remove worst entries when over capacity
    fn prune(&mut self) {
        // Sort by fitness (lower is better) and keep best
        self.entries
            .sort_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap());
        self.entries.truncate(self.max_size);
    }

    /// Get the number of entries
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Get all entries (for clustering)
    pub fn entries(&self) -> &[HallOfFameEntry] {
        &self.entries
    }

    /// Add top-K chromosomes from a population
    pub fn add_from_population(&mut self, chromosomes: &[chromosome::Chromosome], generation: u32) {
        // Get indices sorted by fitness (best first)
        let mut indices: Vec<usize> = (0..chromosomes.len()).collect();
        indices.sort_by(|&a, &b| {
            chromosomes[a]
                .fitness
                .partial_cmp(&chromosomes[b].fitness)
                .unwrap()
        });

        // Try to add top-K
        let mut added = 0;
        for &idx in &indices {
            if added >= HALL_OF_FAME_TOP_K {
                break;
            }
            let chr = &chromosomes[idx];
            if self.try_add(&chr.genes, chr.fitness, generation) {
                added += 1;
            }
        }
    }
}

impl Default for HallOfFame {
    fn default() -> Self {
        Self::new()
    }
}

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
    let mut generation_count = 0u64;
    let mut best_score_history: Vec<f64> = Vec::new();
    let mut generations_without_improvement = 0;
    let mut last_best_score = f64::MAX;

    // Initialize Hall of Fame for tracking diverse solutions
    let mut hall_of_fame = HallOfFame::new();

    println!(
        "{} Starting evolution for {} generations",
        "🧬".bold(),
        MAX_GENERATIONS
    );
    while generation_count < MAX_GENERATIONS {
        pop.eval_fitness();

        // Track diverse solutions in Hall of Fame
        hall_of_fame.add_from_population(&pop.chromosomes, generation_count as u32);

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

    // Finish progress bar if not already finished
    if !progress.is_finished() {
        progress.finish();
    }

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

        // Final evaluation to get best models from final population
        pop.eval_fitness();

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

        // Reconstruct structures from Hall of Fame entries
        let hof_entries = hall_of_fame.entries();
        let hof_structures: Vec<(usize, Molecule, f64)> = hof_entries
            .iter()
            .enumerate()
            .map(|(idx, entry)| {
                let ligand = ligand_clone
                    .clone()
                    .rotate(entry.genes[0], entry.genes[1], entry.genes[2])
                    .displace(entry.genes[3], entry.genes[4], entry.genes[5]);
                let complex = combine_molecules(&receptor_clone, &ligand);
                (idx, complex, entry.fitness)
            })
            .collect();

        // Run FCC clustering
        let structures_only: Vec<Molecule> = hof_structures
            .iter()
            .map(|(_, mol, _)| mol.clone())
            .collect();
        let cluster_config = ClusteringConfig::default();
        let clusters = clustering::cluster_structures(&structures_only, &cluster_config);

        println!(
            "  {} Found {} clusters (min size: {})",
            "✓".green(),
            clusters.len().to_string().cyan(),
            cluster_config.min_cluster_size
        );

        // Select models: cluster centers sorted by energy (best first)
        let mut cluster_centers: Vec<(usize, f64, usize)> = clusters
            .iter()
            .map(|cluster| {
                let fitness = hof_entries[cluster.center_idx].fitness;
                (cluster.center_idx, fitness, cluster.size)
            })
            .collect();

        // Sort by fitness (lowest/best energy first)
        cluster_centers.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // Take up to NUM_OUTPUT_MODELS
        let mut selected: Vec<(usize, usize)> = cluster_centers
            .iter()
            .take(NUM_OUTPUT_MODELS)
            .map(|(idx, _, size)| (*idx, *size))
            .collect();

        // If fewer than NUM_OUTPUT_MODELS clusters, fill with best-scored entries
        if selected.len() < NUM_OUTPUT_MODELS {
            let mut all_entries: Vec<(usize, f64)> = hof_entries
                .iter()
                .enumerate()
                .map(|(idx, entry)| (idx, entry.fitness))
                .collect();
            all_entries.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

            let selected_indices: Vec<usize> = selected.iter().map(|(idx, _)| *idx).collect();

            for (idx, _) in all_entries {
                if selected.len() >= NUM_OUTPUT_MODELS {
                    break;
                }
                if !selected_indices.contains(&idx) {
                    selected.push((idx, 1));
                }
            }
        }

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

        for (model_num, (hof_idx, cluster_size)) in selected.iter().enumerate() {
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

        // Sort Hall of Fame entries by fitness (best first)
        let mut ranked_entries: Vec<(usize, f64)> = hof_entries
            .iter()
            .enumerate()
            .map(|(idx, entry)| (idx, entry.fitness))
            .collect();
        ranked_entries.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        for (rank, (hof_idx, _)) in ranked_entries.iter().take(NUM_OUTPUT_MODELS).enumerate() {
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

    println!("\n{}", "✨ Done!".bold().green());
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_combine_molecules() {
        let receptor_model = read_pdb(&"data/2oob_A.pdb".to_string());
        let ligand_model = read_pdb(&"data/2oob_B.pdb".to_string());

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
    fn test_hall_of_fame_new() {
        let hof = HallOfFame::new();
        assert!(hof.is_empty());
        assert_eq!(hof.len(), 0);
    }

    #[test]
    fn test_hall_of_fame_add_entry() {
        let mut hof = HallOfFame::new();
        let genes = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

        let added = hof.try_add(&genes, -100.0, 1);
        assert!(added);
        assert_eq!(hof.len(), 1);
    }

    #[test]
    fn test_hall_of_fame_rejects_duplicates() {
        let mut hof = HallOfFame::new();
        let genes = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

        hof.try_add(&genes, -100.0, 1);
        let added = hof.try_add(&genes, -100.0, 2);

        assert!(!added, "Should reject identical genes");
        assert_eq!(hof.len(), 1);
    }

    #[test]
    fn test_hall_of_fame_accepts_different_genes() {
        let mut hof = HallOfFame::new();

        // Different rotations (large difference)
        let genes1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let genes2 = [PI, PI, PI, 0.0, 0.0, 0.0];

        hof.try_add(&genes1, -100.0, 1);
        let added = hof.try_add(&genes2, -90.0, 2);

        assert!(added, "Should accept different genes");
        assert_eq!(hof.len(), 2);
    }

    #[test]
    fn test_hall_of_fame_accepts_different_translations() {
        let mut hof = HallOfFame::new();

        // Same rotation, different translation
        let genes1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let genes2 = [0.0, 0.0, 0.0, 10.0, 10.0, 10.0];

        hof.try_add(&genes1, -100.0, 1);
        let added = hof.try_add(&genes2, -90.0, 2);

        assert!(added, "Should accept different translations");
        assert_eq!(hof.len(), 2);
    }

    #[test]
    fn test_hall_of_fame_rejects_similar_genes() {
        let mut hof = HallOfFame::new();

        // Very similar genes (within threshold)
        let genes1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let genes2 = [0.1, 0.1, 0.1, 0.5, 0.5, 0.5]; // Within thresholds

        hof.try_add(&genes1, -100.0, 1);
        let added = hof.try_add(&genes2, -90.0, 2);

        assert!(!added, "Should reject similar genes within threshold");
        assert_eq!(hof.len(), 1);
    }

    #[test]
    fn test_hall_of_fame_prune_keeps_best() {
        let mut hof = HallOfFame::with_capacity(3);

        // Add 5 entries with different genes
        let base_translation = [0.0, 0.0, 0.0];
        for i in 0..5 {
            let genes = [
                (i as f64) * 1.0, // Vary rotation significantly
                0.0,
                0.0,
                base_translation[0] + (i as f64) * 10.0, // Vary translation
                base_translation[1],
                base_translation[2],
            ];
            hof.try_add(&genes, -(100.0 - i as f64 * 10.0), i as u32);
        }

        // Should have pruned to max_size, keeping best fitness
        assert_eq!(hof.len(), 3);

        // Best entries should be kept (lowest fitness)
        let entries = hof.entries();
        assert!(entries[0].fitness <= entries[1].fitness);
        assert!(entries[1].fitness <= entries[2].fitness);
    }

    #[test]
    fn test_angular_difference_same() {
        let diff = HallOfFame::angular_difference(0.0, 0.0);
        assert!((diff - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_angular_difference_wrap_around() {
        // 0 and 2π should be very close
        let diff = HallOfFame::angular_difference(0.0, 2.0 * PI - 0.1);
        assert!(diff < 0.2, "Should handle wrap-around");
    }

    #[test]
    fn test_angular_difference_opposite() {
        // 0 and π should be π apart
        let diff = HallOfFame::angular_difference(0.0, PI);
        assert!((diff - PI).abs() < 0.001);
    }
}
