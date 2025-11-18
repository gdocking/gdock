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
use constants::{MAX_GENERATIONS, POPULATION_SIZE, ENABLE_EARLY_STOPPING, CONVERGENCE_THRESHOLD, CONVERGENCE_WINDOW, DEFAULT_W_AIR, DEFAULT_W_DESOLV, DEFAULT_W_ELEC, DEFAULT_W_VDW};
use structure::read_pdb;


fn score() {
    let mut filenames = Vec::new();

    // --------------------------------------------
    // TODO: Read the filenames from somwehere here
    filenames.push("data/2oob.pdb");
    // --------------------------------------------

    // Loop over the filenames
    for filename in filenames {
        let (mol_a, mol_b) = scoring::read_complex(filename);

        let restraints = restraints::create_restraints(&mol_a, &mol_b);

        let ratio = fitness::satisfaction_ratio(&restraints, &mol_a, &mol_b);

        println!("Ratio of satisfied restraints: {:.3}", ratio);

        let vdw = fitness::vdw_energy(&mol_a, &mol_b);
        let elec = fitness::elec_energy(&mol_a, &mol_b);
        let desolv = fitness::desolv_energy(&mol_a, &mol_b);
        let air = fitness::air_energy(&restraints, &mol_a, &mol_b);

        // Our score function (matching the GA):
        //  1.0*vdw + 1.0*elec + 1.0*desolv + 10.0*air (optimized weights)
        let score = 1.0 * vdw + 1.0 * elec + 1.0 * desolv + 10.0 * air;

        // Print the results
        println!(
            "{} score: {:.3} vdw: {:.3} elec: {:.3} desolv: {:.3} air: {:.3} rest_sat: {:.1}%",
            filename,
            score,
            vdw,
            elec,
            desolv,
            air,
            ratio * 100.0
        );
    }
}

fn run(
    receptor_file: String,
    ligand_file: String,
    restraint_pairs: Vec<(i32, i32)>,
    reference_file: Option<String>,
    w_vdw: f64,
    w_elec: f64,
    w_desolv: f64,
    w_air: f64,
) {
    println!("=== GDock Run Configuration ===");
    println!("Receptor: {}", receptor_file);
    println!("Ligand: {}", ligand_file);
    println!("Restraint pairs: {} (receptor:ligand)", restraint_pairs.len());
    for (i, (rec, lig)) in restraint_pairs.iter().enumerate() {
        println!("  {}. {}:{}", i + 1, rec, lig);
    }
    println!("Reference structure: {}", reference_file.as_ref().unwrap_or(&"None (score-only mode)".to_string()));
    println!("Energy weights: VDW={:.2}, Elec={:.2}, Desolv={:.2}, AIR={:.2}", 
        w_vdw, w_elec, w_desolv, w_air);
    println!();

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
    let restraints = restraints::create_restraints_from_pairs(
        &receptor,
        &ligand,
        &restraint_pairs,
    );
    let num_restraints = restraints.len();
    println!(
        "Created {} restraints from {} residue pairs",
        num_restraints,
        restraint_pairs.len()
    );

    // Clone the original molecule for potential RMSD calculations
    let orig = ligand.clone();

    let ligand = utils::position_ligand(&receptor, ligand);

    // Start the evaluator (only if reference is provided)
    let evaluator = if let Some(ref_file) = &reference_file {
        println!("Loading reference structure for DockQ calculations...");
        let (_, reference_ligand) = scoring::read_complex(ref_file);
        // Use the INPUT receptor (not reference receptor) since receptor doesn't move
        Some(evaluator::Evaluator::new(receptor.clone(), reference_ligand))
    } else {
        println!("No reference structure provided - running in score-only mode");
        None
    };

    // Clone molecules before moving them into population (needed for PDB output later)
    let receptor_clone = receptor.clone();
    let ligand_clone = ligand.clone();

    // Start the GA
    // Create the initial population
    let mut population =
        population::Population::new(Vec::new(), receptor, ligand, orig, restraints, w_vdw, w_elec, w_desolv, w_air);
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
        let metric_vec = if let Some(ref eval) = evaluator {
            Some(population.eval_metrics(eval))
        } else {
            None
        };

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
            println!("\n⚠ Early stopping triggered: No improvement for {} generations", CONVERGENCE_WINDOW);
            println!("Converged at generation {}", generation_count);
            break;
        }

        // Calculate population averages
        let mean_fitness = population.get_mean_fitness();
        let mean_rest: f64 = 100.0 * (population.chromosomes.iter()
            .map(|c| 1.0 - c.restraint_penalty / num_restraints as f64).sum::<f64>() 
            / population.chromosomes.len() as f64);
        
        let best_rest = 100.0 * (1.0 - best_chr.restraint_penalty / num_restraints as f64);
        
        // Calculate improvement metrics
        let improvement_since_last = if generation_count > 0 {
            let prev = best_score_history[generation_count as usize - 1];
            ((prev - best_fitness) / prev.abs()) * 100.0
        } else {
            0.0
        };
        
        let improvement_from_start = if generation_count > 0 {
            ((best_score_history[0] - best_fitness) / best_score_history[0]) * 100.0
        } else {
            0.0
        };
        
        // Print output based on whether we have metrics or not
        if let Some(ref metrics) = metric_vec {
            let mean_dockq: f64 = metrics.iter().map(|m| m.dockq).sum::<f64>() / metrics.len() as f64;
            let best_metrics = &metrics[best_fitness_idx];
            
            println!("Gen #{:04} | Avg: Score={:.1} DockQ={:.3} Rest={:.0}% | Best: Score={:.1} DockQ={:.3} Rest={:.0}% RMSD={:.2} FNAT={:.3} iRMSD={:.2} | Δ={:.2}% ΔTotal={:.1}%",
                generation_count,
                mean_fitness,
                mean_dockq, 
                mean_rest,
                best_fitness, 
                best_metrics.dockq,
                best_rest,
                best_metrics.rmsd,
                best_metrics.fnat,
                best_metrics.irmsd,
                improvement_since_last,
                improvement_from_start
            );
        } else {
            println!("Gen #{:04} | Avg: Score={:.1} Rest={:.0}% | Best: Score={:.1} Rest={:.0}% | Δ={:.2}% ΔTotal={:.1}%",
                generation_count,
                mean_fitness,
                mean_rest,
                best_fitness, 
                best_rest,
                improvement_since_last,
                improvement_from_start
            );
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

    // Save the best models to disk
    println!("\n=== Saving Models ===");

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

        println!("Saved: best_by_score.pdb, best_by_dockq.pdb");
    } else {
        println!("Saved: best_by_score.pdb");
    }
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
    let matches = Command::new("gdock")
        .version("2.0.0")
        .author("pending")
        .about("pending")
        .subcommand(
            Command::new("run")
                .about("Run the GA docking algorithm")
                .arg(clap::Arg::new("receptor")
                    .long("receptor")
                    .short('r')
                    .value_name("FILE")
                    .help("Receptor PDB file")
                    .required(true))
                .arg(clap::Arg::new("ligand")
                    .long("ligand")
                    .short('l')
                    .value_name("FILE")
                    .help("Ligand PDB file")
                    .required(true))
                .arg(clap::Arg::new("restraints")
                    .long("restraints")
                    .value_name("PAIRS")
                    .help("Comma-separated restraint pairs receptor:ligand (e.g., 10:45,15:50,20:55)")
                    .required(true))
                .arg(clap::Arg::new("reference")
                    .long("reference")
                    .value_name("FILE")
                    .help("Reference PDB file for DockQ calculation (optional)"))
                .arg(clap::Arg::new("w_vdw")
                    .long("w_vdw")
                    .value_name("WEIGHT")
                    .help("Weight for VDW energy term (default: 1.0)")
                    .value_parser(clap::value_parser!(f64)))
                .arg(clap::Arg::new("w_elec")
                    .long("w_elec")
                    .value_name("WEIGHT")
                    .help("Weight for electrostatic energy term (default: 0.5)")
                    .value_parser(clap::value_parser!(f64)))
                .arg(clap::Arg::new("w_desolv")
                    .long("w_desolv")
                    .value_name("WEIGHT")
                    .help("Weight for desolvation energy term (default: 0.5)")
                    .value_parser(clap::value_parser!(f64)))
                .arg(clap::Arg::new("w_air")
                    .long("w_air")
                    .value_name("WEIGHT")
                    .help("Weight for AIR restraint energy term (default: 100.0)")
                    .value_parser(clap::value_parser!(f64)))
        )
        .subcommand(Command::new("score").about("Score the results"))
        .get_matches();

    match matches.subcommand() {
        Some(("run", sub_matches)) => {
            let receptor_file = sub_matches.get_one::<String>("receptor").unwrap().clone();
            let ligand_file = sub_matches.get_one::<String>("ligand").unwrap().clone();
            
            // Parse restraint pairs (format: "rec1:lig1,rec2:lig2,...")
            let restraint_pairs: Vec<(i32, i32)> = sub_matches
                .get_one::<String>("restraints")
                .unwrap()
                .split(',')
                .map(|pair| {
                    let parts: Vec<&str> = pair.trim().split(':').collect();
                    if parts.len() != 2 {
                        panic!("Invalid restraint format: '{}'. Expected format: 'receptor:ligand'", pair);
                    }
                    let rec = parts[0].trim().parse::<i32>()
                        .expect(&format!("Invalid receptor residue number: '{}'", parts[0]));
                    let lig = parts[1].trim().parse::<i32>()
                        .expect(&format!("Invalid ligand residue number: '{}'", parts[1]));
                    (rec, lig)
                })
                .collect();
            
            let reference_file = sub_matches.get_one::<String>("reference").cloned();

            let w_vdw = sub_matches.get_one::<f64>("w_vdw").copied().unwrap_or(DEFAULT_W_VDW);
            let w_elec = sub_matches.get_one::<f64>("w_elec").copied().unwrap_or(DEFAULT_W_ELEC);
            let w_desolv = sub_matches.get_one::<f64>("w_desolv").copied().unwrap_or(DEFAULT_W_DESOLV);
            let w_air = sub_matches.get_one::<f64>("w_air").copied().unwrap_or(DEFAULT_W_AIR);
            
            run(
                receptor_file,
                ligand_file,
                restraint_pairs,
                reference_file,
                w_vdw,
                w_elec,
                w_desolv,
                w_air,
            );
        },
        Some(("score", _)) => score(),
        _ => eprintln!("Please specify a valid subcommand: run or score"),
    }
}
