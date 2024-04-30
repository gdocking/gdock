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

use clap::{App, SubCommand};
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
        // let elec = fitness::elec_energy(&mol_a, &mol_b);
        // let desolv = fitness::desolv_energy(&mol_a, &mol_b);
        // let bsa = fitness::bsa_energy(&mol_a, &mol_b);
        // let air = fitness::air_energy(&mol_a, &mol_b);

        // The haddock-score is:
        //  0.01*vdw + 0.1*elec + 1.0*desolv - 0.01*bsa + 0.01*air
        let score = vdw;

        // Print the results
        println!("{} {} {}", filename, score, vdw);
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

    // --------------------------------------------------------------------------------
    // Clone the original molecule for RMSD calculations
    let orig = ligand.clone();

    let ligand = utils::position_ligand(&receptor, ligand);

    // Start the evaluator
    let evaluator = evaluator::Evaluator::new(receptor.clone(), orig.clone());

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

        // -------------------------------------------------------------------------------
        // Do something with the Metrics here
        let metric_vec = population.eval_metrics(&evaluator);

        // Get the minimum rmsd
        let min_rmsd = metric_vec
            .iter()
            .map(|c| c.rmsd)
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        // let min_rmsd = 0.0;

        // --------------------------------------------------------------------------------

        let min_fitness_chromosome = population.get_min_fittest().fitness;
        // let min_fitness_chromosome = 0.0;

        println!(
            "Gen #{:04} min_fitness: {:.3} min_rmsd {:.3}",
            generation_count, min_fitness_chromosome, min_rmsd
        );

        generation_count += 1;
        population = population.evolve(&mut rng);
    }

    println!("Done");
}

fn main() {
    let matches = App::new("gdock")
        .version("2.0.0")
        .author("pending")
        .about("pending")
        .subcommand(SubCommand::with_name("run").about("Run the GA algorithm"))
        .subcommand(SubCommand::with_name("score").about("Score the results"))
        .get_matches();

    match matches.subcommand_name() {
        Some("run") => run(),
        Some("score") => score(),
        _ => eprintln!("Please specify a valid subcommand: run or score"),
    }
}
