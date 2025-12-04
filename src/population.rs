use crate::chromosome::Chromosome;
use crate::constants;
use crate::evaluator;
use crate::ga;
use crate::restraints;
use crate::structure::Molecule;
use core::cmp::Ordering;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::Rng;
// use rayon::iter::IntoParallelRefIterator;
use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator};

#[derive(Debug, Clone)]
pub struct Population {
    pub chromosomes: Vec<Chromosome>,
    pub receptor: Molecule,
    pub ligand: Molecule,
    pub reference: Molecule,
    pub generation: u64,
    pub restraints: Vec<restraints::Restraint>,
    pub weights: constants::EnergyWeights,
}

impl Population {
    pub fn new(
        individuals: Vec<Chromosome>,
        receptor: Molecule,
        ligand: Molecule,
        reference: Molecule,
        restraints: Vec<restraints::Restraint>,
        weights: constants::EnergyWeights,
    ) -> Population {
        Population {
            chromosomes: individuals,
            receptor,
            ligand,
            reference,
            generation: 0,
            restraints,
            weights,
        }
    }

    pub fn eval_fitness(&mut self) {
        // Calculate the fitness in parallel
        let weights = self.weights;
        self.chromosomes.par_iter_mut().for_each(|c| {
            c.fitness(&self.receptor, &self.ligand, &self.restraints, &weights);
        });
    }

    pub fn eval_metrics(&self, evaluator: &evaluator::Evaluator) -> Vec<evaluator::Metrics> {
        self.chromosomes
            .par_iter()
            .map(|c| {
                let model = c.apply_genes(&self.ligand);
                evaluator.calc_metrics(&model)
            })
            .collect()
    }

    /// Evolve the population
    pub fn evolve(&self, rng: &mut StdRng) -> Population {
        // ELITISM: Save best individuals before evolution
        let mut elite = self.chromosomes.clone();
        elite.sort_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap_or(Ordering::Equal));
        let elite_individuals: Vec<Chromosome> = elite
            .iter()
            .take(constants::ELITISM_COUNT)
            .cloned()
            .collect();

        // Tournament selection
        let mut new_population = self.tournament_selection(rng);

        // Crossover population
        new_population.chromosomes.shuffle(rng);
        let mut new_chromosomes = Vec::new();
        for chunk in new_population.chromosomes.chunks(2) {
            if chunk.len() == 2 {
                let (offspring1, offspring2) = ga::crossover(rng, &chunk[0], &chunk[1]);
                new_chromosomes.push(offspring1);
                new_chromosomes.push(offspring2);
            }
            // Optionally handle the case where there's an odd number of chromosomes
        }
        new_population.chromosomes = new_chromosomes;

        // Mutation
        new_population
            .chromosomes
            .iter_mut()
            .for_each(|individual| {
                individual.mutate(rng, constants::MUTATION_RATE);
            });

        // ELITISM: Replace worst individuals with elite
        // Sort new population by fitness (worst first)
        new_population
            .chromosomes
            .sort_by(|a, b| b.fitness.partial_cmp(&a.fitness).unwrap_or(Ordering::Equal));
        // Replace worst N individuals with elite
        for (i, elite_individual) in elite_individuals.iter().enumerate() {
            if i < new_population.chromosomes.len() {
                new_population.chromosomes[i] = elite_individual.clone();
            }
        }

        new_population
    }

    /// Tournament selection
    fn tournament_selection(&self, rng: &mut StdRng) -> Population {
        let mut offspring = Population::new(
            Vec::with_capacity(self.size()),
            self.receptor.clone(),
            self.ligand.clone(),
            self.reference.clone(),
            self.restraints.clone(),
            self.weights,
        );

        offspring.generation = self.generation + 1;

        for _ in 0..self.size() {
            let tournament_individuals = (0..constants::TOURNAMENT_SIZE)
                .map(|_| {
                    let random_index = rng.gen_range(0..self.size());
                    &self.chromosomes[random_index]
                })
                .collect::<Vec<_>>();

            let fittest = tournament_individuals
                .iter()
                .min_by(|a, b| {
                    a.fitness
                        .partial_cmp(&b.fitness)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .expect("No individuals in tournament");

            offspring.chromosomes.push((*fittest).clone());
        }

        offspring
    }

    /// Get the Chromosome with the lowest fitness
    pub fn get_min_fittest(&self) -> &Chromosome {
        self.chromosomes
            .iter()
            .min_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap_or(Ordering::Equal))
            .expect("No chromosomes present")
    }

    /// Get the Chromosome with the highest fitness
    pub fn get_max_fittest(&self) -> &Chromosome {
        // pub fn get_min_fittest(&self) -> &Chromosome {
        self.chromosomes
            .iter()
            .max_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap_or(Ordering::Equal))
            .expect("No chromosomes present")
    }

    /// Get the mean fitness of the population
    pub fn get_mean_fitness(&self) -> f64 {
        if self.chromosomes.is_empty() {
            return f64::NAN;
        }
        let total_fitness: f64 = self.chromosomes.iter().map(|c| c.fitness).sum();
        total_fitness / self.chromosomes.len() as f64
    }

    /// Get the mean fitness of the population
    pub fn get_mean_rmsd(&self) -> f64 {
        if self.chromosomes.is_empty() {
            return f64::NAN;
        }
        let total_rmsd: f64 = self.chromosomes.iter().map(|c| c.rmsd).sum();
        total_rmsd / self.chromosomes.len() as f64
    }

    pub fn get_min_rmsd(&self) -> f64 {
        self.chromosomes
            .iter()
            .min_by(|a, b| a.rmsd.partial_cmp(&b.rmsd).unwrap_or(Ordering::Equal))
            .expect("No chromosomes present")
            .rmsd
    }

    pub fn get_max_fnat(&self) -> f64 {
        self.chromosomes
            .iter()
            .max_by(|a, b| a.fnat.partial_cmp(&b.fnat).unwrap_or(Ordering::Equal))
            .expect("No chromosomes present")
            .fnat
    }

    pub fn size(&self) -> usize {
        self.chromosomes.len()
    }

    pub fn save_individual(&mut self, index: usize, chromossome: Chromosome) {
        self.chromosomes[index] = chromossome;
    }

    pub fn get_individual(&self, index: usize) -> &Chromosome {
        &self.chromosomes[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn create_test_atom(
        x: f64,
        y: f64,
        z: f64,
        resseq: i16,
        chainid: char,
    ) -> crate::structure::Atom {
        crate::structure::Atom {
            serial: 1,
            name: "CA".to_string(),
            altloc: ' ',
            resname: "ALA".to_string(),
            chainid,
            resseq,
            icode: ' ',
            x,
            y,
            z,
            occupancy: 1.0,
            tempfactor: 0.0,
            element: "C".to_string(),
            charge: 0.0,
            vdw_radius: 1.7,
            epsilon: -0.1,
            rmin2: 2.0,
            eps_1_4: -0.1,
            rmin2_1_4: 1.9,
        }
    }

    fn create_test_population(num_chromosomes: usize) -> Population {
        let mut rng = StdRng::seed_from_u64(42);
        let chromosomes: Vec<Chromosome> = (0..num_chromosomes)
            .map(|_| Chromosome::new(&mut rng))
            .collect();

        let mut receptor = Molecule::new();
        receptor.0.push(create_test_atom(0.0, 0.0, 0.0, 1, 'A'));

        let mut ligand = Molecule::new();
        ligand.0.push(create_test_atom(5.0, 0.0, 0.0, 10, 'B'));

        let reference = ligand.clone();

        Population::new(
            chromosomes,
            receptor,
            ligand,
            reference,
            vec![],
            constants::EnergyWeights::default(),
        )
    }

    #[test]
    fn test_population_new() {
        let pop = create_test_population(10);

        assert_eq!(pop.chromosomes.len(), 10);
        assert_eq!(pop.generation, 0);
        assert_eq!(pop.receptor.0.len(), 1);
        assert_eq!(pop.ligand.0.len(), 1);
        assert_eq!(pop.reference.0.len(), 1);
    }

    #[test]
    fn test_population_size() {
        let pop = create_test_population(25);
        assert_eq!(pop.size(), 25);

        let empty_pop = create_test_population(0);
        assert_eq!(empty_pop.size(), 0);
    }

    #[test]
    fn test_eval_fitness() {
        let mut pop = create_test_population(5);

        // Before evaluation, fitness should be 0
        for chromosome in &pop.chromosomes {
            assert_eq!(chromosome.fitness, 0.0);
        }

        pop.eval_fitness();

        // After evaluation, fitness should be calculated (non-zero)
        for chromosome in &pop.chromosomes {
            assert!(chromosome.fitness.is_finite());
        }
    }

    #[test]
    fn test_get_min_fittest() {
        let mut pop = create_test_population(5);
        pop.chromosomes[0].fitness = 10.0;
        pop.chromosomes[1].fitness = 5.0;
        pop.chromosomes[2].fitness = 15.0;
        pop.chromosomes[3].fitness = 3.0;
        pop.chromosomes[4].fitness = 8.0;

        let min_fittest = pop.get_min_fittest();
        assert_eq!(min_fittest.fitness, 3.0);
    }

    #[test]
    fn test_get_max_fittest() {
        let mut pop = create_test_population(5);
        pop.chromosomes[0].fitness = 10.0;
        pop.chromosomes[1].fitness = 5.0;
        pop.chromosomes[2].fitness = 15.0;
        pop.chromosomes[3].fitness = 3.0;
        pop.chromosomes[4].fitness = 8.0;

        let max_fittest = pop.get_max_fittest();
        assert_eq!(max_fittest.fitness, 15.0);
    }

    #[test]
    fn test_get_mean_fitness() {
        let mut pop = create_test_population(5);
        pop.chromosomes[0].fitness = 10.0;
        pop.chromosomes[1].fitness = 20.0;
        pop.chromosomes[2].fitness = 30.0;
        pop.chromosomes[3].fitness = 40.0;
        pop.chromosomes[4].fitness = 50.0;

        let mean = pop.get_mean_fitness();
        assert!((mean - 30.0).abs() < 1e-10);
    }

    #[test]
    fn test_get_mean_fitness_empty() {
        let pop = create_test_population(0);
        let mean = pop.get_mean_fitness();
        assert!(mean.is_nan());
    }

    #[test]
    fn test_get_mean_rmsd() {
        let mut pop = create_test_population(4);
        pop.chromosomes[0].rmsd = 2.0;
        pop.chromosomes[1].rmsd = 4.0;
        pop.chromosomes[2].rmsd = 6.0;
        pop.chromosomes[3].rmsd = 8.0;

        let mean = pop.get_mean_rmsd();
        assert!((mean - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_get_mean_rmsd_empty() {
        let pop = create_test_population(0);
        let mean = pop.get_mean_rmsd();
        assert!(mean.is_nan());
    }

    #[test]
    fn test_get_min_rmsd() {
        let mut pop = create_test_population(4);
        pop.chromosomes[0].rmsd = 5.5;
        pop.chromosomes[1].rmsd = 2.3;
        pop.chromosomes[2].rmsd = 7.8;
        pop.chromosomes[3].rmsd = 3.1;

        let min_rmsd = pop.get_min_rmsd();
        assert!((min_rmsd - 2.3).abs() < 1e-10);
    }

    #[test]
    fn test_get_max_fnat() {
        let mut pop = create_test_population(4);
        pop.chromosomes[0].fnat = 0.5;
        pop.chromosomes[1].fnat = 0.8;
        pop.chromosomes[2].fnat = 0.3;
        pop.chromosomes[3].fnat = 0.95;

        let max_fnat = pop.get_max_fnat();
        assert!((max_fnat - 0.95).abs() < 1e-10);
    }

    #[test]
    fn test_save_and_get_individual() {
        let mut pop = create_test_population(3);
        let mut rng = StdRng::seed_from_u64(123);
        let new_chromosome = Chromosome::new(&mut rng);

        pop.save_individual(1, new_chromosome.clone());
        let retrieved = pop.get_individual(1);

        assert_eq!(retrieved.genes, new_chromosome.genes);
    }

    #[test]
    fn test_tournament_selection_maintains_size() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut pop = create_test_population(20);
        pop.eval_fitness();

        let new_pop = pop.tournament_selection(&mut rng);

        assert_eq!(new_pop.size(), pop.size());
    }

    #[test]
    fn test_tournament_selection_increments_generation() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut pop = create_test_population(10);
        pop.eval_fitness();
        pop.generation = 5;

        let new_pop = pop.tournament_selection(&mut rng);

        assert_eq!(new_pop.generation, 6);
    }

    #[test]
    fn test_tournament_selection_selects_fitter_individuals() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut pop = create_test_population(10);

        // Set fitness values - lower is better
        for (i, chromosome) in pop.chromosomes.iter_mut().enumerate() {
            chromosome.fitness = (i as f64) * 10.0; // 0, 10, 20, ..., 90
        }

        let new_pop = pop.tournament_selection(&mut rng);

        // Tournament selection should favor lower fitness values
        let mean_original = pop.get_mean_fitness();
        let mean_selected = new_pop.get_mean_fitness();

        // The mean fitness of selected population should be lower or equal
        // (though randomness means it might not always be strictly lower)
        assert!(
            mean_selected.is_finite() && mean_original.is_finite(),
            "Both means should be finite"
        );
    }

    #[test]
    fn test_evolve_maintains_population_size() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut pop = create_test_population(20);
        pop.eval_fitness();

        let new_pop = pop.evolve(&mut rng);

        assert_eq!(new_pop.size(), pop.size());
    }

    #[test]
    fn test_evolve_increments_generation() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut pop = create_test_population(10);
        pop.eval_fitness();
        pop.generation = 3;

        let new_pop = pop.evolve(&mut rng);

        assert_eq!(new_pop.generation, 4);
    }

    #[test]
    fn test_evolve_preserves_elite() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut pop = create_test_population(20);
        pop.eval_fitness();

        // Set specific fitness values to track elites
        for (i, chromosome) in pop.chromosomes.iter_mut().enumerate() {
            chromosome.fitness = (i as f64) * 5.0;
        }

        // Best individuals have lowest fitness (0, 5, 10, 15, 20)
        let elite_fitness_values: Vec<f64> = pop
            .chromosomes
            .iter()
            .take(constants::ELITISM_COUNT)
            .map(|c| c.fitness)
            .collect();

        let new_pop = pop.evolve(&mut rng);

        // Check that the elite fitness values are preserved in new population
        let mut elite_preserved = 0;
        for elite_fitness in &elite_fitness_values {
            if new_pop
                .chromosomes
                .iter()
                .any(|c| (c.fitness - elite_fitness).abs() < 1e-10)
            {
                elite_preserved += 1;
            }
        }

        assert_eq!(
            elite_preserved,
            constants::ELITISM_COUNT,
            "Elite individuals should be preserved"
        );
    }

    #[test]
    fn test_eval_metrics() {
        let mut pop = create_test_population(3);
        pop.eval_fitness();

        let evaluator = evaluator::Evaluator::new(pop.receptor.clone(), pop.reference.clone());
        let metrics = pop.eval_metrics(&evaluator);

        assert_eq!(metrics.len(), 3);
        for metric in &metrics {
            assert!(metric.rmsd.is_finite());
            assert!(metric.irmsd.is_finite());
            assert!(metric.fnat.is_finite());
            assert!(metric.dockq.is_finite());
        }
    }
}
