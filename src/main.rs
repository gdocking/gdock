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
use constants::{MAX_GENERATIONS, POPULATION_SIZE};
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

fn run() {
    // --------------------------------------------------------------------------------
    // Define number of threads either with the function below or with the
    // environment variable RAYON_NUM_THREADS
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(1)
    //     .build_global()
    //     .unwrap();
    // --------------------------------------------------------------------------------

    let input_receptor = &"data/A.pdb".to_string();
    let input_ligand = &"data/B.pdb".to_string();

    let receptor = read_pdb(input_receptor);
    let ligand = read_pdb(input_ligand);

    // --------------------------------------------------------------------------------
    // DEVELOPMENT //
    // --------------------------------------------------------------------------------
    let restraints = restraints::create_restraints(&receptor, &ligand);
    let num_restraints = restraints.len();
    println!(
        "Created {} restraints from native structure",
        num_restraints
    );

    // --------------------------------------------------------------------------------
    // Clone the original molecule for RMSD calculations
    let orig = ligand.clone();

    let ligand = utils::position_ligand(&receptor, ligand);

    // Start the evaluator
    let evaluator = evaluator::Evaluator::new(receptor.clone(), orig.clone());

    // Clone molecules before moving them into population (needed for PDB output later)
    let receptor_clone = receptor.clone();
    let ligand_clone = ligand.clone();
    let orig_clone = orig.clone();

    // Start the GA
    // Create the initial population
    let mut population =
        population::Population::new(Vec::new(), receptor, ligand, orig, restraints);
    let mut rng = StdRng::seed_from_u64(constants::RANDOM_SEED);
    for _ in 0..POPULATION_SIZE {
        let c = chromosome::Chromosome::new(&mut rng);
        population.chromosomes.push(c);
    }

    // Evolve the population
    let mut generation_count = 0;
    while generation_count < MAX_GENERATIONS {
        population.eval_fitness();

        // Calculate metrics for all chromosomes
        let metric_vec = population.eval_metrics(&evaluator);

        // Find the best fitness chromosome and its corresponding RMSD
        let best_fitness_idx = population
            .chromosomes
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.fitness.partial_cmp(&b.fitness).unwrap())
            .map(|(idx, _)| idx)
            .unwrap();

        let best_fitness = population.chromosomes[best_fitness_idx].fitness;
        let best_fitness_rmsd = metric_vec[best_fitness_idx].rmsd;
        let best_fitness_fnat = metric_vec[best_fitness_idx].fnat;
        let best_fitness_dockq = metric_vec[best_fitness_idx].dockq;

        // Find the best DockQ and its corresponding fitness
        let best_dockq_idx = metric_vec
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.dockq.partial_cmp(&b.dockq).unwrap())
            .map(|(idx, _)| idx)
            .unwrap();

        let best_dockq = metric_vec[best_dockq_idx].dockq;
        let best_dockq_fitness = population.chromosomes[best_dockq_idx].fitness;

        // Calculate population statistics
        let mean_fitness = population.get_mean_fitness();
        let mean_rmsd: f64 =
            metric_vec.iter().map(|m| m.rmsd).sum::<f64>() / metric_vec.len() as f64;

        // Get energy components for best fitness chromosome
        let best_chr = &population.chromosomes[best_fitness_idx];
        let best_dockq_chr = &population.chromosomes[best_dockq_idx];

        // Calculate restraint satisfaction as percentage
        let best_score_restraints_sat =
            100.0 * (1.0 - best_chr.restraint_penalty / num_restraints as f64);
        let best_dockq_restraints_sat =
            100.0 * (1.0 - best_dockq_chr.restraint_penalty / num_restraints as f64);

        // Print detailed progress every 10 generations
        if generation_count % 10 == 0 {
            println!("Gen #{:04} | BestScore: {:.1} (vdw={:.1}, elec={:.1}, desolv={:.1}, air={:.1}) RMSD={:.2}, Rest={:.0}%, FNAT={:.3}, DockQ={:.3}",
                generation_count, best_fitness, best_chr.vdw, best_chr.elec, best_chr.desolv, best_chr.air,
                best_fitness_rmsd, best_score_restraints_sat, best_fitness_fnat, best_fitness_dockq);
            println!("        | BestDockQ: {:.3} (Score={:.1}, RMSD={:.2}, Rest={:.0}%) | PopAvg: score={:.1}, rmsd={:.2}",
                best_dockq, best_dockq_fitness, metric_vec[best_dockq_idx].rmsd, best_dockq_restraints_sat, mean_fitness, mean_rmsd);
        } else {
            println!("Gen #{:04} | Score: {:.1} [RMSD {:.2}, Rest {:.0}%] | DockQ: {:.3} [Score {:.1}, RMSD {:.2}, Rest {:.0}%]",
                generation_count, best_fitness, best_fitness_rmsd, best_score_restraints_sat,
                best_dockq, best_dockq_fitness, metric_vec[best_dockq_idx].rmsd, best_dockq_restraints_sat);
        }

        generation_count += 1;
        population = population.evolve(&mut rng);
    }

    // Final summary
    println!("\n=== Final Results ===");
    population.eval_fitness();
    let final_metrics = population.eval_metrics(&evaluator);

    let best_fitness_idx = population
        .chromosomes
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.fitness.partial_cmp(&b.fitness).unwrap())
        .map(|(idx, _)| idx)
        .unwrap();

    let best_dockq_idx = final_metrics
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.dockq.partial_cmp(&b.dockq).unwrap())
        .map(|(idx, _)| idx)
        .unwrap();

    let final_best_score = &population.chromosomes[best_fitness_idx];
    let final_best_dockq = &population.chromosomes[best_dockq_idx];
    let final_best_score_sat =
        100.0 * (1.0 - final_best_score.restraint_penalty / num_restraints as f64);
    let final_best_dockq_sat =
        100.0 * (1.0 - final_best_dockq.restraint_penalty / num_restraints as f64);

    println!("Best by Score: fitness={:.2} (vdw={:.1}, elec={:.1}, desolv={:.1}, air={:.1}), RMSD={:.2}, Restraints={:.0}%, FNAT={:.3}, DockQ={:.3}, Rank={}",
        final_best_score.fitness,
        final_best_score.vdw,
        final_best_score.elec,
        final_best_score.desolv,
        final_best_score.air,
        final_metrics[best_fitness_idx].rmsd,
        final_best_score_sat,
        final_metrics[best_fitness_idx].fnat,
        final_metrics[best_fitness_idx].dockq,
        final_metrics[best_fitness_idx].rank());

    println!("Best by DockQ: fitness={:.2} (vdw={:.1}, elec={:.1}, desolv={:.1}, air={:.1}), RMSD={:.2}, Restraints={:.0}%, FNAT={:.3}, DockQ={:.3}, Rank={}",
        final_best_dockq.fitness,
        final_best_dockq.vdw,
        final_best_dockq.elec,
        final_best_dockq.desolv,
        final_best_dockq.air,
        final_metrics[best_dockq_idx].rmsd,
        final_best_dockq_sat,
        final_metrics[best_dockq_idx].fnat,
        final_metrics[best_dockq_idx].dockq,
        final_metrics[best_dockq_idx].rank());

    // Save the best models to disk for visual inspection
    println!("\n=== Saving Models to Disk ===");

    // Apply transformations and save best-by-score model
    let best_score_ligand = final_best_score.apply_genes(&ligand_clone);
    let best_score_complex = combine_molecules(&receptor_clone, &best_score_ligand);
    let output_best_score = "best_by_score.pdb".to_string();
    structure::write_pdb(&best_score_complex, &output_best_score);
    println!("Saved best by score to: {}", output_best_score);
    println!(
        "  Score={:.2}, RMSD={:.2}, Rest={:.0}%, FNAT={:.3}, DockQ={:.3}",
        final_best_score.fitness,
        final_metrics[best_fitness_idx].rmsd,
        final_best_score_sat,
        final_metrics[best_fitness_idx].fnat,
        final_metrics[best_fitness_idx].dockq
    );

    // Apply transformations and save best-by-DockQ model
    let best_dockq_ligand = final_best_dockq.apply_genes(&ligand_clone);
    let best_dockq_complex = combine_molecules(&receptor_clone, &best_dockq_ligand);
    let output_best_dockq = "best_by_dockq.pdb".to_string();
    structure::write_pdb(&best_dockq_complex, &output_best_dockq);
    println!("Saved best by DockQ to: {}", output_best_dockq);
    println!(
        "  Score={:.2}, RMSD={:.2}, Rest={:.0}%, FNAT={:.3}, DockQ={:.3}",
        final_best_dockq.fitness,
        final_metrics[best_dockq_idx].rmsd,
        final_best_dockq_sat,
        final_metrics[best_dockq_idx].fnat,
        final_metrics[best_dockq_idx].dockq
    );

    // Also save the native structure for comparison
    let native_complex = combine_molecules(&receptor_clone, &orig_clone);
    let output_native = "native_reference.pdb".to_string();
    structure::write_pdb(&native_complex, &output_native);
    println!("Saved native reference to: {}", output_native);

    println!("\nDone - You can now visualize the structures in PyMOL/ChimeraX");
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
        .subcommand(Command::new("run").about("Run the GA algorithm"))
        .subcommand(Command::new("score").about("Score the results"))
        .get_matches();

    match matches.subcommand() {
        Some(("run", _)) => run(),
        Some(("score", _)) => score(),
        _ => eprintln!("Please specify a valid subcommand: run or score"),
    }
}
