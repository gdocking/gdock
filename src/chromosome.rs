use super::constants;
use super::constants::MAX_DISPLACEMENT;
use crate::fitness;
use crate::restraints;
use crate::structure;
use core::f64::consts::PI;
use rand::rngs::StdRng;

use rand::Rng;

#[derive(Debug, Clone)]
pub struct Chromosome {
    pub genes: Vec<f64>,
    pub fitness: f64,
    // Energy components for debugging
    pub vdw: f64,
    pub elec: f64,
    pub desolv: f64,
    pub air: f64,
    pub restraint_penalty: f64,
    // Metrics for analysis
    pub rmsd: f64,
    pub fnat: f64,
}

impl Chromosome {
    pub fn new(rng: &mut StdRng) -> Chromosome {
        let alpha = rng.gen_range(0.0..=2.0 * PI);
        let beta = rng.gen_range(0.0..=2.0 * PI);
        let gamma = rng.gen_range(0.0..=2.0 * PI);

        let displacement_x = rng.gen_range(-MAX_DISPLACEMENT..=MAX_DISPLACEMENT);
        let displacement_y = rng.gen_range(-MAX_DISPLACEMENT..=MAX_DISPLACEMENT);
        let displacement_z = rng.gen_range(-MAX_DISPLACEMENT..=MAX_DISPLACEMENT);

        Chromosome {
            genes: vec![
                alpha,
                beta,
                gamma,
                displacement_x,
                displacement_y,
                displacement_z,
            ],
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

    pub fn mutate(&mut self, rng: &mut StdRng, mutation_chance: f64) {
        use rand_distr::{Distribution, Normal};

        // Creep mutation: add small perturbations instead of random replacement
        // This (should) allow local refinement around good solutions
        for i in 0..self.genes.len() {
            if rng.gen_range(0.0..1.0) <= mutation_chance {
                if i < 3 {
                    // Rotations (alpha, beta, gamma): perturbation ±10° (0.174 radians)
                    // Allows fine-tuning of orientation
                    let perturbation = Normal::new(0.0, 0.174).unwrap().sample(rng);
                    self.genes[i] = (self.genes[i] + perturbation).rem_euclid(2.0 * PI);
                } else {
                    // Translations (x, y, z): perturbation ±1.0 Å
                    // Allows fine-tuning of position
                    let perturbation = Normal::new(0.0, 1.0).unwrap().sample(rng);
                    self.genes[i] =
                        (self.genes[i] + perturbation).clamp(-MAX_DISPLACEMENT, MAX_DISPLACEMENT);
                }
            }
        }
    }

    pub fn apply_genes(&self, ligand: &structure::Molecule) -> structure::Molecule {
        let lig = ligand.clone();
        lig.rotate(self.genes[0], self.genes[1], self.genes[2])
            .displace(self.genes[3], self.genes[4], self.genes[5])
    }

    pub fn fitness(
        &mut self,
        receptor: &structure::Molecule,
        ligand: &structure::Molecule,
        restraints: &[restraints::Restraint],
        weights: &constants::EnergyWeights,
        evaluator: Option<&crate::evaluator::Evaluator>,
    ) -> f64 {
        let target_ligand = self.apply_genes(ligand);

        // Debug mode: use ONLY negative DockQ as fitness (lower is better, so we negate)
        // This validates that the sampling strategy can find the native structure
        // WITHOUT using restraints or energy terms - pure quality metric optimization
        if let Some(eval) = evaluator {
            let metrics = eval.calc_metrics(&target_ligand);
            self.fitness = -metrics.dockq; // Negate because GA minimizes fitness
                                           // Set energy components to zero since they're not used in debug mode
            self.vdw = 0.0;
            self.elec = 0.0;
            self.desolv = 0.0;
            self.air = 0.0;
            self.restraint_penalty = 0.0;
            return self.fitness;
        }

        // Normal mode: calculate all energy components
        // vdw 1.0 - prevent clashes and guide packing
        // elec 0.5 - moderate electrostatic influence
        // desolv 0.5 - burial
        // air 100.0 -> let restraints guide id
        self.vdw = fitness::vdw_energy(receptor, &target_ligand);
        self.elec = fitness::elec_energy(receptor, &target_ligand);
        self.desolv = fitness::desolv_energy(receptor, &target_ligand);
        self.air = fitness::air_energy(restraints, receptor, &target_ligand);

        // Calculate restraint satisfaction for monitoring
        let restraints_ratio = fitness::satisfaction_ratio(restraints, receptor, &target_ligand);
        self.restraint_penalty = (1.0 - restraints_ratio) * restraints.len() as f64;

        // Information-driven docking score with configurable weights
        let score = weights.vdw * self.vdw
            + weights.elec * self.elec
            + weights.desolv * self.desolv
            + weights.air * self.air;

        self.fitness = score;
        self.fitness
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use std::f64::consts::PI;

    fn create_test_atom(name: &str, x: f64, y: f64, z: f64) -> structure::Atom {
        structure::Atom {
            serial: 1,
            name: name.to_string(),
            altloc: ' ',
            resname: "ALA".to_string(),
            chainid: 'A',
            resseq: 1,
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

    #[test]
    fn test_chromosome_new() {
        let mut rng = StdRng::seed_from_u64(42);
        let chromosome = Chromosome::new(&mut rng);

        // Should have 6 genes (3 rotations + 3 translations)
        assert_eq!(chromosome.genes.len(), 6);

        // Rotation angles should be in range [0, 2π]
        for i in 0..3 {
            assert!(chromosome.genes[i] >= 0.0);
            assert!(chromosome.genes[i] <= 2.0 * PI);
        }

        // Translations should be in range [-MAX_DISPLACEMENT, MAX_DISPLACEMENT]
        for i in 3..6 {
            assert!(chromosome.genes[i] >= -MAX_DISPLACEMENT);
            assert!(chromosome.genes[i] <= MAX_DISPLACEMENT);
        }

        // Initial fitness values should be zero
        assert_eq!(chromosome.fitness, 0.0);
        assert_eq!(chromosome.vdw, 0.0);
        assert_eq!(chromosome.elec, 0.0);
    }

    #[test]
    fn test_chromosome_mutate_changes_genes() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut chromosome = Chromosome::new(&mut rng);
        let original_genes = chromosome.genes.clone();

        // Mutate with high rate to ensure some changes
        chromosome.mutate(&mut rng, 1.0); // 100% mutation rate

        // At least some genes should have changed
        let changed_count = chromosome
            .genes
            .iter()
            .zip(original_genes.iter())
            .filter(|(new_gene, old_gene)| (*new_gene - *old_gene).abs() > 1e-10)
            .count();

        assert!(
            changed_count > 0,
            "Mutation should change at least one gene"
        );
    }

    #[test]
    fn test_chromosome_mutate_keeps_bounds() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut chromosome = Chromosome::new(&mut rng);

        // Mutate many times
        for _ in 0..100 {
            chromosome.mutate(&mut rng, 0.5);

            // Check rotation bounds [0, 2π]
            for i in 0..3 {
                assert!(
                    chromosome.genes[i] >= 0.0 && chromosome.genes[i] <= 2.0 * PI,
                    "Rotation gene {} out of bounds: {}",
                    i,
                    chromosome.genes[i]
                );
            }

            // Check translation bounds [-MAX_DISPLACEMENT, MAX_DISPLACEMENT]
            for i in 3..6 {
                assert!(
                    chromosome.genes[i] >= -MAX_DISPLACEMENT
                        && chromosome.genes[i] <= MAX_DISPLACEMENT,
                    "Translation gene {} out of bounds: {}",
                    i,
                    chromosome.genes[i]
                );
            }
        }
    }

    #[test]
    fn test_chromosome_mutate_with_zero_rate() {
        let mut rng = StdRng::seed_from_u64(42);
        let mut chromosome = Chromosome::new(&mut rng);
        let original_genes = chromosome.genes.clone();

        chromosome.mutate(&mut rng, 0.0); // 0% mutation rate

        // Genes should be unchanged
        for (new_gene, old_gene) in chromosome.genes.iter().zip(original_genes.iter()) {
            assert_eq!(new_gene, old_gene);
        }
    }

    #[test]
    fn test_apply_genes_translation_only() {
        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 0.0, 0.0, 0.0));

        let chromosome = Chromosome {
            genes: vec![0.0, 0.0, 0.0, 5.0, 10.0, 15.0], // No rotation, just translation
            fitness: 0.0,
            vdw: 0.0,
            elec: 0.0,
            desolv: 0.0,
            air: 0.0,
            restraint_penalty: 0.0,
            rmsd: 0.0,
            fnat: 0.0,
        };

        let transformed = chromosome.apply_genes(&ligand);

        // Position should be translated by (5, 10, 15)
        assert!((transformed.0[0].x - 5.0).abs() < 1e-6);
        assert!((transformed.0[0].y - 10.0).abs() < 1e-6);
        assert!((transformed.0[0].z - 15.0).abs() < 1e-6);
    }

    #[test]
    fn test_apply_genes_preserves_molecule_size() {
        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 0.0, 0.0, 0.0));
        ligand.0.push(create_test_atom("CB", 1.0, 0.0, 0.0));
        ligand.0.push(create_test_atom("C", 2.0, 0.0, 0.0));

        let mut rng = StdRng::seed_from_u64(42);
        let chromosome = Chromosome::new(&mut rng);

        let transformed = chromosome.apply_genes(&ligand);

        assert_eq!(transformed.0.len(), ligand.0.len());
    }

    #[test]
    fn test_fitness_calculation() {
        // Create simple test molecules
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom("CA", 0.0, 0.0, 0.0));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 5.0, 0.0, 0.0));

        let restraints = vec![]; // No restraints
        let weights = constants::EnergyWeights::default();

        let mut chromosome = Chromosome {
            genes: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // Identity transformation
            fitness: 0.0,
            vdw: 0.0,
            elec: 0.0,
            desolv: 0.0,
            air: 0.0,
            restraint_penalty: 0.0,
            rmsd: 0.0,
            fnat: 0.0,
        };

        let fitness = chromosome.fitness(&receptor, &ligand, &restraints, &weights, None);

        // Fitness should be calculated (non-zero for this setup)
        assert!(fitness.is_finite(), "Fitness should be a finite number");
        assert_eq!(
            chromosome.fitness, fitness,
            "Fitness should be stored in chromosome"
        );
    }

    #[test]
    fn test_fitness_stores_energy_components() {
        let mut receptor = structure::Molecule::new();
        receptor.0.push(create_test_atom("CA", 0.0, 0.0, 0.0));

        let mut ligand = structure::Molecule::new();
        ligand.0.push(create_test_atom("CA", 3.0, 0.0, 0.0));

        let restraints = vec![];
        let weights = constants::EnergyWeights::default();

        let mut chromosome = Chromosome {
            genes: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            fitness: 0.0,
            vdw: 0.0,
            elec: 0.0,
            desolv: 0.0,
            air: 0.0,
            restraint_penalty: 0.0,
            rmsd: 0.0,
            fnat: 0.0,
        };

        chromosome.fitness(&receptor, &ligand, &restraints, &weights, None);

        // Energy components should be calculated and stored
        assert!(
            chromosome.vdw.is_finite(),
            "VDW energy should be calculated"
        );
        assert!(
            chromosome.elec.is_finite(),
            "Electrostatic energy should be calculated"
        );
        assert!(
            chromosome.desolv.is_finite(),
            "Desolvation energy should be calculated"
        );
    }
}
