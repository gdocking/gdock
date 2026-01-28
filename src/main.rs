pub mod chromosome;
pub mod clustering;
pub mod commands;
pub mod constants;
pub mod evaluator;
pub mod fitness;
pub mod ga;
pub mod population;
pub mod restraints;
pub mod scoring;
pub mod structure;
pub mod toppar;
pub mod utils;

use clap::Command;
use constants::{DEFAULT_W_AIR, DEFAULT_W_DESOLV, DEFAULT_W_ELEC, DEFAULT_W_VDW};
use std::fs::File;
use std::io::prelude::*;

/// Helper to parse restraint pairs from string format "rec1:lig1,rec2:lig2,..."
fn parse_restraints(restraints_str: &str) -> Result<Vec<(i32, i32)>, String> {
    restraints_str
        .split(',')
        .filter(|pair| !pair.trim().is_empty())
        .map(|pair| {
            let parts: Vec<&str> = pair.trim().split(':').collect();
            if parts.len() != 2 {
                return Err(format!(
                    "Invalid restraint format: '{}'. Expected format: 'receptor:ligand'",
                    pair
                ));
            }
            let rec = parts[0]
                .trim()
                .parse::<i32>()
                .map_err(|_| format!("Invalid receptor residue number: '{}'", parts[0]))?;
            let lig = parts[1]
                .trim()
                .parse::<i32>()
                .map_err(|_| format!("Invalid ligand residue number: '{}'", parts[1]))?;
            Ok((rec, lig))
        })
        .collect()
}

/// Helper functions to parse restraints pair from a file
fn parse_restraints_file(restraints_file_path: &str) -> Result<Vec<(i32, i32)>, String> {
    let mut file = File::open(restraints_file_path).map_err(|e| {
        format!(
            "Cannot open restraints file '{}': {}",
            restraints_file_path, e
        )
    })?;
    let mut contents = String::new();
    file.read_to_string(&mut contents).map_err(|e| {
        format!(
            "Cannot read restraints file '{}': {}",
            restraints_file_path, e
        )
    })?;
    parse_restraints(&contents)
}

fn main() {
    const VERSION: &str = env!("CARGO_PKG_VERSION");

    // Common arguments for receptor and ligand
    let receptor_arg = clap::Arg::new("receptor")
        .long("receptor")
        .short('r')
        .value_name("FILE")
        .help("Receptor PDB file")
        .required(true);

    let ligand_arg = clap::Arg::new("ligand")
        .long("ligand")
        .short('l')
        .value_name("FILE")
        .help("Ligand PDB file")
        .required(true);

    // Weight arguments (shared by run and score)
    let weight_args = [
        clap::Arg::new("w_vdw")
            .long("w_vdw")
            .value_name("WEIGHT")
            .help("Weight for VDW energy term")
            .value_parser(clap::value_parser!(f64)),
        clap::Arg::new("w_elec")
            .long("w_elec")
            .value_name("WEIGHT")
            .help("Weight for electrostatic energy term")
            .value_parser(clap::value_parser!(f64)),
        clap::Arg::new("w_desolv")
            .long("w_desolv")
            .value_name("WEIGHT")
            .help("Weight for desolvation energy term")
            .value_parser(clap::value_parser!(f64)),
        clap::Arg::new("w_air")
            .long("w_air")
            .value_name("WEIGHT")
            .help("Weight for AIR restraint energy term")
            .value_parser(clap::value_parser!(f64)),
    ];

    let matches = Command::new("gdock")
        .version(VERSION)
        .author("Rodrigo V. Honorato")
        .about("Fast information-driven protein-protein docking using genetic algorithms")
        .subcommand_required(true)
        .subcommand(
            Command::new("run")
                .about("Run the genetic algorithm docking")
                .arg(receptor_arg.clone())
                .arg(ligand_arg.clone())
                .arg(
                    clap::Arg::new("restraints")
                        .long("restraints")
                        .value_name("PAIRS")
                        .help("Comma-separated restraint pairs receptor:ligand (e.g., 10:45,15:50)")
                        .required(true),
                )
                .arg(
                    clap::Arg::new("reference")
                        .long("reference")
                        .value_name("FILE")
                        .help("Reference PDB file for DockQ calculation"),
                )
                .arg(
                    clap::Arg::new("debug")
                        .long("debug")
                        .help("Debug mode: use DockQ as fitness (requires --reference)")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    clap::Arg::new("output-dir")
                        .long("output-dir")
                        .short('o')
                        .value_name("DIR")
                        .help("Output directory for results (default: current directory)"),
                )
                .arg(
                    clap::Arg::new("no-clust")
                        .long("no-clust")
                        .help("Disable clustering, output best_by_score and best_by_dockq only")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    clap::Arg::new("nproc")
                        .long("nproc")
                        .short('n')
                        .value_name("NUM")
                        .help("Number of processors to use (default: total - 2)")
                        .value_parser(clap::value_parser!(usize)),
                )
                .args(weight_args.clone()),
        )
        .subcommand(
            Command::new("score")
                .about("Score structures without running the GA")
                .arg(receptor_arg.clone())
                .arg(ligand_arg.clone())
                .arg(
                    clap::Arg::new("restraints")
                        .long("restraints")
                        .value_name("PAIRS")
                        .help("Comma-separated restraint pairs receptor:ligand (optional)"),
                )
                .arg(
                    clap::Arg::new("reference")
                        .long("reference")
                        .value_name("FILE")
                        .help("Reference PDB file for DockQ calculation"),
                )
                .args(weight_args),
        )
        .subcommand(
            Command::new("restraints")
                .about("Generate restraints from interface contacts")
                .arg(receptor_arg)
                .arg(ligand_arg)
                .arg(
                    clap::Arg::new("cutoff")
                        .long("cutoff")
                        .value_name("ANGSTROMS")
                        .help("Distance cutoff for interface detection (default: 5.0)")
                        .value_parser(clap::value_parser!(f64)),
                ),
        )
        .get_matches();

    match matches.subcommand() {
        Some(("run", sub_m)) => {
            // Configure thread pool based on --nproc
            let total_cpus = std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1);
            let default_threads = if total_cpus > 2 { total_cpus - 2 } else { 1 };
            let requested_threads = sub_m.get_one::<usize>("nproc").copied();
            let num_threads = match requested_threads {
                Some(n) if n > total_cpus => {
                    eprintln!(
                        "Warning: requested {} threads but only {} available, using {}",
                        n, total_cpus, default_threads
                    );
                    default_threads
                }
                Some(0) => {
                    eprintln!(
                        "Warning: --nproc must be at least 1, using {}",
                        default_threads
                    );
                    default_threads
                }
                Some(n) => n,
                None => default_threads,
            };
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .unwrap_or_else(|e| eprintln!("Warning: could not set thread pool: {}", e));

            let receptor_file = sub_m.get_one::<String>("receptor").unwrap().clone();
            let ligand_file = sub_m.get_one::<String>("ligand").unwrap().clone();
            let reference_file = sub_m.get_one::<String>("reference").cloned();
            let debug_mode = sub_m.get_flag("debug");

            // Validate debug mode
            if debug_mode && reference_file.is_none() {
                eprintln!("Error: --debug requires --reference to be specified");
                std::process::exit(1);
            }

            let restraints_arg = sub_m.get_one::<String>("restraints").unwrap();
            let restraint_pairs = if std::path::Path::new(restraints_arg).is_file() {
                parse_restraints_file(restraints_arg)
            } else {
                parse_restraints(restraints_arg)
            }
            .unwrap_or_else(|e| {
                eprintln!("Error: {}", e);
                std::process::exit(1);
            });

            let weights = constants::EnergyWeights::new(
                sub_m
                    .get_one::<f64>("w_vdw")
                    .copied()
                    .unwrap_or(DEFAULT_W_VDW),
                sub_m
                    .get_one::<f64>("w_elec")
                    .copied()
                    .unwrap_or(DEFAULT_W_ELEC),
                sub_m
                    .get_one::<f64>("w_desolv")
                    .copied()
                    .unwrap_or(DEFAULT_W_DESOLV),
                sub_m
                    .get_one::<f64>("w_air")
                    .copied()
                    .unwrap_or(DEFAULT_W_AIR),
            );

            let output_dir = sub_m.get_one::<String>("output-dir").cloned();
            let no_clustering = sub_m.get_flag("no-clust");

            commands::run::run(commands::run::RunConfig {
                receptor_file,
                ligand_file,
                restraint_pairs,
                reference_file,
                weights,
                debug_mode,
                output_dir,
                no_clustering,
            });
        }
        Some(("score", sub_m)) => {
            let receptor_file = sub_m.get_one::<String>("receptor").unwrap().clone();
            let ligand_file = sub_m.get_one::<String>("ligand").unwrap().clone();
            let reference_file = sub_m.get_one::<String>("reference").cloned();

            let restraint_pairs = sub_m.get_one::<String>("restraints").map(|s| {
                let result = if std::path::Path::new(s).is_file() {
                    parse_restraints_file(s)
                } else {
                    parse_restraints(s)
                };
                result.unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                })
            });

            let weights = constants::EnergyWeights::new(
                sub_m
                    .get_one::<f64>("w_vdw")
                    .copied()
                    .unwrap_or(DEFAULT_W_VDW),
                sub_m
                    .get_one::<f64>("w_elec")
                    .copied()
                    .unwrap_or(DEFAULT_W_ELEC),
                sub_m
                    .get_one::<f64>("w_desolv")
                    .copied()
                    .unwrap_or(DEFAULT_W_DESOLV),
                sub_m
                    .get_one::<f64>("w_air")
                    .copied()
                    .unwrap_or(DEFAULT_W_AIR),
            );

            commands::score::score(
                receptor_file,
                ligand_file,
                restraint_pairs,
                reference_file,
                weights,
            );
        }
        Some(("restraints", sub_m)) => {
            let receptor_file = sub_m.get_one::<String>("receptor").unwrap().clone();
            let ligand_file = sub_m.get_one::<String>("ligand").unwrap().clone();
            let cutoff = sub_m.get_one::<f64>("cutoff").copied().unwrap_or(5.0);

            commands::restraints::generate_restraints(receptor_file, ligand_file, cutoff);
        }
        _ => unreachable!("subcommand_required prevents this"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_restraints_single() {
        let pairs = parse_restraints("10:45").unwrap();
        assert_eq!(pairs, vec![(10, 45)]);
    }

    #[test]
    fn test_parse_restraints_multiple() {
        let pairs = parse_restraints("10:45,15:50,20:55").unwrap();
        assert_eq!(pairs, vec![(10, 45), (15, 50), (20, 55)]);
    }

    #[test]
    fn test_parse_restraints_with_spaces() {
        let pairs = parse_restraints("10:45, 15:50 , 20:55").unwrap();
        assert_eq!(pairs, vec![(10, 45), (15, 50), (20, 55)]);
    }

    #[test]
    fn test_parse_restraints_invalid_format() {
        let result = parse_restraints("10-45");
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Invalid restraint format"));
    }

    #[test]
    fn test_parse_restraints_invalid_number() {
        let result = parse_restraints("abc:45");
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .contains("Invalid receptor residue number"));
    }
}
