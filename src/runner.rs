//! Top-level GA runner — the single source of truth for the evolution loop.
use std::collections::HashSet;

use rand::rngs::StdRng;

use crate::clustering::{self, ClusteringConfig};
use crate::constants::{
    CONVERGENCE_THRESHOLD, CONVERGENCE_WINDOW, ENABLE_EARLY_STOPPING, NUM_OUTPUT_MODELS,
};
use crate::hall_of_fame::{HallOfFame, HallOfFameEntry};
use crate::population::Population;
use crate::structure::{combine_molecules, Molecule};

/// Result returned by [`run_ga`] after the evolution loop completes.
pub struct GaResult {
    pub hall_of_fame: HallOfFame,
    /// Number of the last generation that was evaluated (0-indexed).
    pub generations_run: u64,
    /// `true` if the loop terminated because of the early-stopping criterion.
    pub converged_early: bool,
    /// The population after the final `eval_fitness` call.
    pub final_population: Population,
}

/// Run the genetic algorithm evolution loop.
///
/// `on_generation` is called once per generation **after** fitness evaluation and
/// Hall-of-Fame update, but **before** the early-stopping check. This lets callers
/// (e.g. the CLI) display per-generation progress without duplicating the core loop
/// logic. Pass `|_, _| {}` when no reporting is needed.
pub fn run_ga<F>(
    mut pop: Population,
    rng: &mut StdRng,
    max_generations: u64,
    mut on_generation: F,
) -> GaResult
where
    F: FnMut(u64, &Population),
{
    let mut hall_of_fame = HallOfFame::new();
    let mut generation_count = 0u64;
    let mut generations_without_improvement = 0u64;
    let mut last_best_score = f64::MAX;
    let mut converged_early = false;

    while generation_count < max_generations {
        pop.eval_fitness();
        hall_of_fame.add_from_population(&pop.chromosomes);

        let best_fitness = pop.get_min_fittest().fitness;

        if generation_count > 0 {
            let improvement = (last_best_score - best_fitness) / last_best_score.abs();
            if improvement < CONVERGENCE_THRESHOLD {
                generations_without_improvement += 1;
            } else {
                generations_without_improvement = 0;
            }
        }
        last_best_score = best_fitness;

        on_generation(generation_count, &pop);

        if ENABLE_EARLY_STOPPING && generations_without_improvement >= CONVERGENCE_WINDOW {
            converged_early = true;
            break;
        }

        generation_count += 1;
        pop = pop.evolve(rng);
    }

    // Final evaluation so callers always receive a fully-evaluated population.
    pop.eval_fitness();
    hall_of_fame.add_from_population(&pop.chromosomes);

    GaResult {
        hall_of_fame,
        generations_run: generation_count,
        converged_early,
        final_population: pop,
    }
}

/// Output of [`select_models`]: the HoF indices to use for each output set.
pub struct SelectedModels {
    /// `(hof_idx, cluster_size)` — cluster-centre models sorted by fitness (best first).
    pub clustered: Vec<(usize, usize)>,
    /// `hof_idx` — top entries sorted purely by fitness (best first).
    pub ranked: Vec<usize>,
}

/// Select the output models from a Hall-of-Fame by FCC clustering.
///
/// Returns up to `NUM_OUTPUT_MODELS` cluster-centre entries (diverse poses) and
/// `NUM_OUTPUT_MODELS` entries ranked purely by fitness. When there are fewer
/// clusters than `NUM_OUTPUT_MODELS`, the gap is filled with the next-best
/// unselected HoF entries.
pub fn select_models(
    hof_entries: &[HallOfFameEntry],
    receptor: &Molecule,
    ligand: &Molecule,
) -> SelectedModels {
    // Build receptor+ligand complexes for every HoF entry.
    let complexes: Vec<Molecule> = hof_entries
        .iter()
        .map(|e| {
            let docked = ligand
                .clone()
                .rotate(e.genes[0], e.genes[1], e.genes[2])
                .displace(e.genes[3], e.genes[4], e.genes[5]);
            combine_molecules(receptor, &docked)
        })
        .collect();

    // FCC clustering.
    let cluster_config = ClusteringConfig::default();
    let clusters = clustering::cluster_structures(&complexes, &cluster_config);

    // Pick cluster centres sorted by fitness.
    let mut cluster_centers: Vec<(usize, f64, usize)> = clusters
        .iter()
        .map(|c| (c.center_idx, hof_entries[c.center_idx].fitness, c.size))
        .collect();
    cluster_centers.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut clustered: Vec<(usize, usize)> = cluster_centers
        .iter()
        .take(NUM_OUTPUT_MODELS)
        .map(|(idx, _, size)| (*idx, *size))
        .collect();

    // Fill to NUM_OUTPUT_MODELS with the best unselected HoF entries.
    let mut used: HashSet<usize> = clustered.iter().map(|(i, _)| *i).collect();
    if clustered.len() < NUM_OUTPUT_MODELS {
        let mut remaining: Vec<(usize, f64)> = hof_entries
            .iter()
            .enumerate()
            .filter(|(i, _)| !used.contains(i))
            .map(|(i, e)| (i, e.fitness))
            .collect();
        remaining.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        for (idx, _) in remaining.iter().take(NUM_OUTPUT_MODELS - clustered.len()) {
            clustered.push((*idx, 1));
            used.insert(*idx);
        }
    }

    // Ranked purely by fitness.
    let mut ranked_by_fitness: Vec<(usize, f64)> = hof_entries
        .iter()
        .enumerate()
        .map(|(i, e)| (i, e.fitness))
        .collect();
    ranked_by_fitness
        .sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let ranked: Vec<usize> = ranked_by_fitness
        .iter()
        .take(NUM_OUTPUT_MODELS)
        .map(|(idx, _)| *idx)
        .collect();

    SelectedModels { clustered, ranked }
}
