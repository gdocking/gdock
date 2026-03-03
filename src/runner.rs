//! Top-level GA runner — the single source of truth for the evolution loop.
use rand::rngs::StdRng;

use crate::constants::{CONVERGENCE_THRESHOLD, CONVERGENCE_WINDOW, ENABLE_EARLY_STOPPING};
use crate::hall_of_fame::HallOfFame;
use crate::population::Population;

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
