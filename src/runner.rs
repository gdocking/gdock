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
    hof_capacity: usize,
    mut on_generation: F,
) -> GaResult
where
    F: FnMut(u64, &Population),
{
    let mut hall_of_fame = HallOfFame::with_capacity(hof_capacity);
    let mut generation_count = 0u64;
    let mut generations_without_improvement = 0u64;
    let mut last_best_score = f64::MAX;
    let mut converged_early = false;

    while generation_count < max_generations {
        pop.eval_fitness();
        hall_of_fame.add_from_population(&pop.chromosomes);

        let best_fitness = pop.get_min_fittest().fitness;

        if generation_count > 0 {
            let improvement = if last_best_score.abs() < f64::EPSILON {
                0.0
            } else {
                (last_best_score - best_fitness) / last_best_score.abs()
            };
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
    // Skipped on early-stop: the last loop iteration already evaluated and added to HoF.
    if !converged_early {
        pop.eval_fitness();
        hall_of_fame.add_from_population(&pop.chromosomes);
    }

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
    ranked_by_fitness.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let ranked: Vec<usize> = ranked_by_fitness
        .iter()
        .take(NUM_OUTPUT_MODELS)
        .map(|(idx, _)| *idx)
        .collect();

    SelectedModels { clustered, ranked }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chromosome::Chromosome;
    use crate::constants::{EnergyWeights, HALL_OF_FAME_MAX_SIZE};
    use crate::population::Population;
    use crate::structure::Molecule;
    use rand::SeedableRng;

    fn minimal_population(n: usize) -> Population {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let chromosomes: Vec<Chromosome> = (0..n).map(|_| Chromosome::new(&mut rng)).collect();
        let mut receptor = Molecule::new();
        receptor.0.push(crate::structure::Atom {
            serial: 1,
            name: "CA".to_string(),
            altloc: ' ',
            resname: "ALA".to_string(),
            chainid: 'A',
            resseq: 1,
            icode: ' ',
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            tempfactor: 0.0,
            element: "C".to_string(),
            charge: 0.0,
            vdw_radius: 1.7,
            epsilon: -0.1,
            rmin2: 2.0,
            eps_1_4: -0.1,
            rmin2_1_4: 1.9,
        });
        let mut ligand = Molecule::new();
        ligand.0.push(crate::structure::Atom {
            serial: 2,
            name: "CA".to_string(),
            altloc: ' ',
            resname: "ALA".to_string(),
            chainid: 'B',
            resseq: 1,
            icode: ' ',
            x: 5.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            tempfactor: 0.0,
            element: "C".to_string(),
            charge: 0.0,
            vdw_radius: 1.7,
            epsilon: -0.1,
            rmin2: 2.0,
            eps_1_4: -0.1,
            rmin2_1_4: 1.9,
        });
        let reference = ligand.clone();
        Population::new(
            chromosomes,
            receptor,
            ligand,
            reference,
            vec![],
            EnergyWeights::default(),
            None,
        )
    }

    #[test]
    fn test_run_ga_hof_capacity_respected() {
        let pop = minimal_population(10);
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let capacity = 2;
        let result = run_ga(pop, &mut rng, 5, capacity, |_, _| {});
        assert!(
            result.hall_of_fame.len() <= capacity,
            "HoF should not exceed hof_capacity"
        );
    }

    #[test]
    fn test_run_ga_hof_never_exceeds_capacity() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let result = run_ga(minimal_population(10), &mut rng, 10, 1, |_, _| {});
        assert!(
            result.hall_of_fame.len() <= 1,
            "HoF should not exceed capacity=1"
        );

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let result = run_ga(
            minimal_population(10),
            &mut rng,
            10,
            HALL_OF_FAME_MAX_SIZE,
            |_, _| {},
        );
        assert!(
            result.hall_of_fame.len() <= HALL_OF_FAME_MAX_SIZE,
            "HoF should not exceed HALL_OF_FAME_MAX_SIZE"
        );
    }
}
