use rand::rngs::StdRng;
use rand::SeedableRng;
use std::cmp::Ordering;
use std::fs;
use std::path::{Path, PathBuf};

use crate::chromosome;
use crate::constants::{
    self, EnergyWeights, HALL_OF_FAME_MAX_SIZE, MAX_GENERATIONS, POPULATION_SIZE,
};
use crate::evaluator;
use crate::hall_of_fame::{HallOfFame, HallOfFameEntry};
use crate::population;
use crate::reporting;
use crate::restraints;
use crate::runner::{run_ga, select_models, GaResult};
use crate::scoring;
use crate::structure::{self, read_pdb, Molecule};
use crate::utils;

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
    pub epoch: bool,
}

pub use crate::structure::combine_molecules;

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
        epoch,
    } = config;
    const VERSION: &str = env!("CARGO_PKG_VERSION");

    let receptor_model = read_pdb(&receptor_file);
    let ligand_model = read_pdb(&ligand_file);
    let receptor = receptor_model.0[0].clone();
    let ligand = ligand_model.0[0].clone();

    let restraints_list =
        restraints::create_ambiguous_restraints_from_pairs(&receptor, &ligand, &restraint_pairs);
    let num_restraints = restraints_list.len() / 2;

    reporting::header(
        VERSION,
        debug_mode,
        &receptor_file,
        &ligand_file,
        &reference_file,
        restraint_pairs.len(),
        num_restraints,
        &weights,
    );

    let orig = ligand.clone();
    let ligand = utils::position_ligand(&receptor, ligand);

    let eval = if let Some(ref_file) = &reference_file {
        let (_, reference_ligand) = scoring::read_complex(ref_file);
        Some(evaluator::Evaluator::new(
            receptor.clone(),
            reference_ligand,
        ))
    } else {
        None
    };

    // debug mode fits against DockQ instead of the energy score
    let debug_evaluator = if debug_mode { eval.clone() } else { None };

    // kept for PDB output after receptor/ligand are moved into Population::new
    let receptor_clone = receptor.clone();
    let ligand_clone = ligand.clone();

    let out_dir = match &output_dir {
        Some(dir) => {
            let path = PathBuf::from(dir);
            fs::create_dir_all(&path).expect("Failed to create output directory");
            path
        }
        None => PathBuf::from("."),
    };

    let epoch_writer = if epoch {
        Some(reporting::EpochWriter::create(&out_dir, eval.is_some()))
    } else {
        None
    };

    let mut rng = StdRng::seed_from_u64(constants::RANDOM_SEED);
    let mut pop = population::Population::new(
        Vec::new(),
        receptor,
        ligand,
        orig,
        restraints_list,
        weights,
        debug_evaluator,
    );
    for _ in 0..POPULATION_SIZE {
        pop.chromosomes.push(chromosome::Chromosome::new(&mut rng));
    }

    let hof_capacity = sampling.unwrap_or(HALL_OF_FAME_MAX_SIZE);
    let ga_result = evolve(
        pop,
        &mut rng,
        eval.as_ref(),
        debug_mode,
        hof_capacity,
        epoch_writer,
    );

    let ctx = OutputContext {
        receptor: &receptor_clone,
        ligand: &ligand_clone,
        eval: eval.as_ref(),
        out_dir: &out_dir,
        debug_mode,
        epoch,
        generations_run: ga_result.generations_run,
        converged_early: ga_result.converged_early,
    };

    if no_clustering {
        write_best_models_output(&ctx, &ga_result.final_population);
    } else {
        write_clustered_output(&ctx, &ga_result.hall_of_fame);
    }

    if sampling.is_some() {
        write_sampling_output(
            &ga_result.hall_of_fame,
            &out_dir,
            &receptor_clone,
            &ligand_clone,
        );
    }

    reporting::done();
}

/// Run the GA loop, reporting live progress and (when enabled) streaming the
/// per-generation trace to `epoch_writer`.
fn evolve(
    pop: population::Population,
    rng: &mut StdRng,
    eval: Option<&evaluator::Evaluator>,
    debug_mode: bool,
    hof_capacity: usize,
    mut epoch_writer: Option<reporting::EpochWriter>,
) -> GaResult {
    let progress = reporting::new_progress_bar(MAX_GENERATIONS);
    reporting::evolution_start(MAX_GENERATIONS, POPULATION_SIZE);

    let ga_result = run_ga(pop, rng, MAX_GENERATIONS, hof_capacity, |gen, pop| {
        let metric_vec = eval.map(|e| pop.eval_metrics(e));

        if let Some(w) = epoch_writer.as_mut() {
            for (i, chr) in pop.chromosomes.iter().enumerate() {
                w.write_row(gen, i, chr, metric_vec.as_ref().map(|m| &m[i]));
            }
        }

        let best_idx = pop
            .chromosomes
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.fitness.partial_cmp(&b.fitness).unwrap())
            .map(|(idx, _)| idx)
            .unwrap();
        let best_fitness = pop.chromosomes[best_idx].fitness;

        progress.set_position(gen);

        let message = match metric_vec.as_ref() {
            Some(metrics) => {
                reporting::progress_message(debug_mode, best_fitness, Some(metrics), best_idx, 0.0)
            }
            None => reporting::progress_message(
                debug_mode,
                best_fitness,
                None,
                best_idx,
                pop.get_mean_fitness(),
            ),
        };
        progress.set_message(message);
    });

    if let Some(w) = epoch_writer {
        w.finish();
    }

    reporting::finish_progress(
        &progress,
        ga_result.generations_run,
        ga_result.converged_early,
    );

    ga_result
}

/// Shared inputs threaded through the output writers, borrowed from `run()`.
struct OutputContext<'a> {
    receptor: &'a Molecule,
    ligand: &'a Molecule,
    eval: Option<&'a evaluator::Evaluator>,
    out_dir: &'a Path,
    debug_mode: bool,
    epoch: bool,
    generations_run: u64,
    converged_early: bool,
}

/// `--no-clust`: write `best_by_score.pdb`, and (with a reference)
/// `best_by_dockq.pdb` plus `metrics.tsv`.
fn write_best_models_output(ctx: &OutputContext, pop: &population::Population) {
    let best_fitness_idx = pop
        .chromosomes
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.fitness.partial_cmp(&b.fitness).unwrap())
        .map(|(idx, _)| idx)
        .unwrap();

    let final_best_score = &pop.chromosomes[best_fitness_idx];

    reporting::saving_results();

    let best_score_ligand = final_best_score.apply_genes(ctx.ligand);
    let best_score_complex = combine_molecules(ctx.receptor, &best_score_ligand);
    let best_score_path = ctx.out_dir.join("best_by_score.pdb");
    structure::write_pdb(
        &best_score_complex,
        best_score_path.to_string_lossy().as_ref(),
    );

    let Some(e) = ctx.eval else {
        reporting::saved_file(&best_score_path);
        if ctx.epoch {
            reporting::saved_file(&ctx.out_dir.join("epoch.tsv.gz"));
        }
        return;
    };

    let final_metrics = pop.eval_metrics(e);
    let best_dockq_idx = final_metrics
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.dockq.partial_cmp(&b.dockq).unwrap())
        .map(|(idx, _)| idx)
        .unwrap();

    let final_best_dockq = &pop.chromosomes[best_dockq_idx];
    let best_dockq_ligand = final_best_dockq.apply_genes(ctx.ligand);
    let best_dockq_complex = combine_molecules(ctx.receptor, &best_dockq_ligand);
    let best_dockq_path = ctx.out_dir.join("best_by_dockq.pdb");
    structure::write_pdb(
        &best_dockq_complex,
        best_dockq_path.to_string_lossy().as_ref(),
    );

    let best_score_metrics = &final_metrics[best_fitness_idx];
    let best_dockq_metrics = &final_metrics[best_dockq_idx];

    reporting::final_metrics(best_score_metrics, best_dockq_metrics);

    let metrics_path = ctx.out_dir.join("metrics.tsv");
    let mut metrics_file = fs::File::create(&metrics_path).expect("Failed to create metrics file");
    reporting::write_noclust_metrics_header(&mut metrics_file);
    reporting::write_noclust_metrics_row(
        &mut metrics_file,
        "best_by_score",
        best_score_metrics,
        final_best_score.fitness,
        ctx.generations_run,
        ctx.converged_early,
    );
    reporting::write_noclust_metrics_row(
        &mut metrics_file,
        "best_by_dockq",
        best_dockq_metrics,
        final_best_dockq.fitness,
        ctx.generations_run,
        ctx.converged_early,
    );

    reporting::saved_file(&best_score_path);
    reporting::saved_file(&best_dockq_path);
    reporting::saved_file(&metrics_path);
    if ctx.epoch {
        reporting::saved_file(&ctx.out_dir.join("epoch.tsv.gz"));
    }
}

/// Default path: cluster the Hall of Fame and write the diverse `model_*` set
/// plus the top-scoring `ranked_*` set, all into one `metrics.tsv`.
fn write_clustered_output(ctx: &OutputContext, hall_of_fame: &HallOfFame) {
    reporting::clustering(hall_of_fame.len());

    let hof_entries = hall_of_fame.entries();
    let selected = select_models(hof_entries, ctx.receptor, ctx.ligand);

    let metrics_path = ctx.out_dir.join("metrics.tsv");
    let mut metrics_file = fs::File::create(&metrics_path).expect("Failed to create metrics file");
    reporting::write_clustered_metrics_header(&mut metrics_file, ctx.eval.is_some());

    reporting::model_section_header(
        "📊 Clustered models",
        "· diverse poses",
        !ctx.debug_mode,
        ctx.eval.is_some(),
        true,
    );
    for (model_num, (hof_idx, cluster_size)) in selected.clustered.iter().enumerate() {
        let model_name = format!("model_{}", model_num + 1);
        write_model(
            ctx,
            &mut metrics_file,
            &hof_entries[*hof_idx],
            &model_name,
            Some(*cluster_size),
        );
    }

    reporting::model_section_header(
        "📊 Ranked models",
        "· best score",
        !ctx.debug_mode,
        ctx.eval.is_some(),
        false,
    );
    for (rank, hof_idx) in selected.ranked.iter().enumerate() {
        let model_name = format!("ranked_{}", rank + 1);
        write_model(
            ctx,
            &mut metrics_file,
            &hof_entries[*hof_idx],
            &model_name,
            None,
        );
    }

    let n_models = selected.clustered.len() + selected.ranked.len();
    reporting::summary(n_models, ctx.epoch, ctx.out_dir);
}

/// Write one model's PDB, its `metrics.tsv` row, and its on-screen table row.
/// `cluster_size` is `Some` for clustered models, `None` for ranked ones.
fn write_model(
    ctx: &OutputContext,
    metrics_file: &mut fs::File,
    entry: &HallOfFameEntry,
    model_name: &str,
    cluster_size: Option<usize>,
) {
    let ligand = ctx
        .ligand
        .clone()
        .rotate(entry.genes[0], entry.genes[1], entry.genes[2])
        .displace(entry.genes[3], entry.genes[4], entry.genes[5]);
    let complex = combine_molecules(ctx.receptor, &ligand);

    let pdb_path = ctx.out_dir.join(format!("{model_name}.pdb"));
    structure::write_pdb(&complex, pdb_path.to_string_lossy().as_ref());

    let metrics = ctx.eval.map(|e| e.calc_metrics(&ligand));
    reporting::write_clustered_metrics_row(
        metrics_file,
        model_name,
        cluster_size,
        entry.fitness,
        metrics.as_ref(),
        ctx.generations_run,
        ctx.converged_early,
    );

    let score = (!ctx.debug_mode).then_some(entry.fitness);
    let dockq = metrics.as_ref().map(|m| m.dockq);
    reporting::model_row(model_name, score, dockq, cluster_size);
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
    reporting::write_sampling_header(&mut tsv);

    reporting::sampling_header();

    for (rank, entry) in sorted.iter().enumerate() {
        let model_name = format!("gdock_{}", rank + 1);
        let docked_ligand = ligand
            .clone()
            .rotate(entry.genes[0], entry.genes[1], entry.genes[2])
            .displace(entry.genes[3], entry.genes[4], entry.genes[5]);
        let complex = combine_molecules(receptor, &docked_ligand);
        let pdb_path = sampling_dir.join(format!("{}.pdb", model_name));
        structure::write_pdb(&complex, pdb_path.to_string_lossy().as_ref());
        reporting::write_sampling_row(&mut tsv, &model_name, entry);
    }

    reporting::sampling_done(sorted.len(), &sampling_dir, &tsv_path);
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
        hof.try_add(
            &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            -100.0,
            1.0,
            2.0,
            3.0,
            4.0,
            0.0,
        );
        hof.try_add(
            &[PI, PI, PI, 10.0, 10.0, 10.0],
            -50.0,
            1.0,
            2.0,
            3.0,
            4.0,
            0.0,
        );

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
        hof.try_add(
            &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            -50.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        );
        hof.try_add(
            &[PI, PI, PI, 10.0, 10.0, 10.0],
            -100.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        );

        let tmp = crate::utils::get_unique_tempdir();
        write_sampling_output(&hof, tmp.as_path(), &receptor, &ligand);

        let tsv = std::fs::read_to_string(tmp.as_path().join("sampling/sampling.tsv")).unwrap();
        let mut lines = tsv.lines();
        assert_eq!(
            lines.next().unwrap(),
            "model\tscore\tvdw\telec\tdesolv\tair\tclash"
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
        hof.try_add(
            &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            -99.0,
            1.1,
            2.2,
            3.3,
            4.4,
            5.5,
        );

        let tmp = crate::utils::get_unique_tempdir();
        write_sampling_output(&hof, tmp.as_path(), &receptor, &ligand);

        let tsv = std::fs::read_to_string(tmp.as_path().join("sampling/sampling.tsv")).unwrap();
        let mut lines = tsv.lines();
        assert_eq!(
            lines.next().unwrap(),
            "model\tscore\tvdw\telec\tdesolv\tair\tclash"
        );

        let data_row = lines.next().unwrap();
        let cols: Vec<&str> = data_row.split('\t').collect();
        assert_eq!(cols.len(), 7, "data rows should have 7 columns");
        assert_eq!(cols[0], "gdock_1");
        assert!((cols[1].parse::<f64>().unwrap() - (-99.0)).abs() < 0.001);
        assert!((cols[2].parse::<f64>().unwrap() - 1.1).abs() < 0.001);
        assert!((cols[3].parse::<f64>().unwrap() - 2.2).abs() < 0.001);
        assert!((cols[4].parse::<f64>().unwrap() - 3.3).abs() < 0.001);
        assert!((cols[5].parse::<f64>().unwrap() - 4.4).abs() < 0.001);
        assert!((cols[6].parse::<f64>().unwrap() - 5.5).abs() < 0.001);
    }
}
