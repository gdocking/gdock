/// This module contains the genetic algorithm implementation.
use super::constants::CROSSOVER_RATE;
use crate::chromosome::Chromosome;
use crate::StdRng;
use rand::Rng;

pub fn crossover(
    rng: &mut StdRng,
    individual1: &Chromosome,
    individual2: &Chromosome,
) -> (Chromosome, Chromosome) {
    let mut new_genes1 = Vec::new();
    let mut new_genes2 = Vec::new();

    // Zip the genes of the individuals together
    let zipped_genes = individual1.genes.iter().zip(individual2.genes.iter());

    // Loop over the zipped genes
    for (gene1, gene2) in zipped_genes {
        // Check if the genes should be swapped
        if rng.gen_range(0.0..1.0) <= CROSSOVER_RATE {
            // Swap the genes
            new_genes1.push(*gene2);
            new_genes2.push(*gene1);
        } else {
            // Don't swap the genes
            new_genes1.push(*gene1);
            new_genes2.push(*gene2);
        }
    }

    let mut new_individual1 = individual1.clone();
    let mut new_individual2 = individual2.clone();

    new_individual1.genes = new_genes1;
    new_individual2.genes = new_genes2;

    (new_individual1, new_individual2)
}
