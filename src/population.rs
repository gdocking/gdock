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
}

impl Population {
    pub fn new(
        individuals: Vec<Chromosome>,
        receptor: Molecule,
        ligand: Molecule,
        reference: Molecule,
        restraints: Vec<restraints::Restraint>,
    ) -> Population {
        Population {
            chromosomes: individuals,
            receptor,
            ligand,
            reference,
            generation: 0,
            restraints,
        }
    }

    pub fn eval_fitness(&mut self) {
        // Calculate the fitness in parallel
        self.chromosomes.par_iter_mut().for_each(|c| {
            c.fitness(&self.receptor, &self.ligand, &self.restraints);
        });
    }

    pub fn eval_metrics(&self, evaluator: &evaluator::Evaluator) -> Vec<evaluator::Metrics> {
        self.chromosomes
            .par_iter()
            .map(|c| {
                let model = c.apply_genes(&self.ligand);
                // let model = self.ligand.clone();
                evaluator.calc_metrics(&model)
            })
            .collect()
    }

    /// Evolve the population
    pub fn evolve(&self, rng: &mut StdRng) -> Population {
        // ELITISM: Save best individuals before evolution
        let mut elite = self.chromosomes.clone();
        elite.sort_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap_or(Ordering::Equal));
        let elite_individuals: Vec<Chromosome> = elite.iter().take(constants::ELITISM_COUNT).cloned().collect();
        
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
        new_population.chromosomes.sort_by(|a, b| b.fitness.partial_cmp(&a.fitness).unwrap_or(Ordering::Equal));
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
            return std::f64::NAN;
        }
        let total_fitness: f64 = self.chromosomes.iter().map(|c| c.fitness).sum();
        total_fitness / self.chromosomes.len() as f64
    }

    /// Get the mean fitness of the population
    pub fn get_mean_rmsd(&self) -> f64 {
        if self.chromosomes.is_empty() {
            return std::f64::NAN;
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
