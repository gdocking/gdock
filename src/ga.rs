/// This module contains the genetic algorithm implementation.
use super::constants::CROSSOVER_RATE;
use crate::chromosome::Chromosome;
use rand::rngs::StdRng;
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn create_test_chromosome(genes: Vec<f64>) -> Chromosome {
        Chromosome {
            genes,
            fitness: 0.0,
            vdw: 0.0,
            elec: 0.0,
            desolv: 0.0,
            air: 0.0,
            restraint_penalty: 0.0,
            rmsd: 0.0,
            fnat: 0.0,
        }
    }

    #[test]
    fn test_crossover_produces_two_offspring() {
        let mut rng = StdRng::seed_from_u64(42);
        let parent1 = create_test_chromosome(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let parent2 = create_test_chromosome(vec![7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);

        let (offspring1, offspring2) = crossover(&mut rng, &parent1, &parent2);

        assert_eq!(offspring1.genes.len(), 6);
        assert_eq!(offspring2.genes.len(), 6);
    }

    #[test]
    fn test_crossover_preserves_gene_count() {
        let mut rng = StdRng::seed_from_u64(42);
        let parent1 = create_test_chromosome(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let parent2 = create_test_chromosome(vec![7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);

        let (offspring1, offspring2) = crossover(&mut rng, &parent1, &parent2);

        assert_eq!(offspring1.genes.len(), parent1.genes.len());
        assert_eq!(offspring2.genes.len(), parent2.genes.len());
    }

    #[test]
    fn test_crossover_genes_come_from_parents() {
        let mut rng = StdRng::seed_from_u64(42);
        let parent1 = create_test_chromosome(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let parent2 = create_test_chromosome(vec![7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);

        let (offspring1, offspring2) = crossover(&mut rng, &parent1, &parent2);

        // Each gene in offspring should come from one of the parents
        for i in 0..6 {
            let gene1 = offspring1.genes[i];
            let gene2 = offspring2.genes[i];

            // Gene must be from either parent1 or parent2
            let from_p1 = gene1 == parent1.genes[i] || gene1 == parent2.genes[i];
            let from_p2 = gene2 == parent1.genes[i] || gene2 == parent2.genes[i];

            assert!(from_p1, "Offspring1 gene {} not from parents", i);
            assert!(from_p2, "Offspring2 gene {} not from parents", i);
        }
    }

    #[test]
    fn test_crossover_swaps_genes() {
        let mut rng = StdRng::seed_from_u64(42);
        let parent1 = create_test_chromosome(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let parent2 = create_test_chromosome(vec![7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);

        let (offspring1, offspring2) = crossover(&mut rng, &parent1, &parent2);

        // When genes are swapped, offspring1 gets gene from parent2 and vice versa
        // Count how many genes were swapped
        let mut swapped_count = 0;
        for i in 0..6 {
            if offspring1.genes[i] == parent2.genes[i] && offspring2.genes[i] == parent1.genes[i] {
                swapped_count += 1;
            }
        }

        // With CROSSOVER_RATE = 0.6, we expect some genes to be swapped
        // With 6 genes and rate 0.6, we expect roughly 3-4 swaps
        assert!(
            swapped_count > 0,
            "Expected at least one gene swap with CROSSOVER_RATE = {}",
            CROSSOVER_RATE
        );
    }

    #[test]
    fn test_crossover_offspring_differ_from_parents() {
        let mut rng = StdRng::seed_from_u64(42);
        let parent1 = create_test_chromosome(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let parent2 = create_test_chromosome(vec![7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);

        let (offspring1, offspring2) = crossover(&mut rng, &parent1, &parent2);

        // With CROSSOVER_RATE = 0.6, offspring should differ from parents
        let offspring1_same = offspring1.genes == parent1.genes;
        let offspring2_same = offspring2.genes == parent2.genes;

        // It's very unlikely both offspring are identical to their parents with rate 0.6
        assert!(
            !offspring1_same || !offspring2_same,
            "At least one offspring should differ from its parent"
        );
    }

    #[test]
    fn test_crossover_complementary_swap() {
        let mut rng = StdRng::seed_from_u64(42);
        let parent1 = create_test_chromosome(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let parent2 = create_test_chromosome(vec![7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);

        let (offspring1, offspring2) = crossover(&mut rng, &parent1, &parent2);

        // For each position, if offspring1 got gene from parent2, then offspring2 got gene from parent1
        for i in 0..6 {
            if offspring1.genes[i] == parent2.genes[i] {
                assert_eq!(
                    offspring2.genes[i], parent1.genes[i],
                    "Complementary swap failed at position {}",
                    i
                );
            } else if offspring1.genes[i] == parent1.genes[i] {
                assert_eq!(
                    offspring2.genes[i], parent2.genes[i],
                    "Complementary swap failed at position {}",
                    i
                );
            }
        }
    }

    #[test]
    fn test_crossover_preserves_other_fields() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut parent1 = create_test_chromosome(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        parent1.fitness = 100.0;
        parent1.vdw = 10.0;

        let parent2 = create_test_chromosome(vec![7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);

        let (offspring1, offspring2) = crossover(&mut rng, &parent1, &parent2);

        // Fitness values should be preserved from parents (will be recalculated later)
        assert_eq!(offspring1.fitness, parent1.fitness);
        assert_eq!(offspring2.fitness, parent2.fitness);
    }
}
