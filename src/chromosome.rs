use super::constants::MAX_DISPLACEMENT;
use crate::fitness;
use crate::restraints;
use crate::structure;
use crate::StdRng;
use core::f64::consts::PI;

use rand::Rng;
// use rayon::prelude::*;

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
        let random_chromossome: Chromosome = Chromosome::new(rng);
        self.genes
            .iter_mut()
            .zip(random_chromossome.genes.iter())
            .for_each(|(gene, random_gene)| {
                if rng.gen_range(0.0..1.0) <= mutation_chance {
                    *gene = *random_gene;
                }
            });
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
    ) -> f64 {
        let target_ligand = self.apply_genes(ligand);

        self.vdw = fitness::vdw_energy(receptor, &target_ligand);
        self.elec = fitness::elec_energy(receptor, &target_ligand);
        self.desolv = fitness::desolv_energy(receptor, &target_ligand);
        self.air = fitness::air_energy(restraints, receptor, &target_ligand);

        // Calculate restraint satisfaction for monitoring
        let restraints_ratio = fitness::satisfaction_ratio(restraints, receptor, &target_ligand);
        self.restraint_penalty = (1.0 - restraints_ratio) * restraints.len() as f64;

        // Information-driven docking score (HADDOCK-like):
        // Physics-based terms guide realistic interactions
        // AIR term guides toward native-like contacts
        // Weights heavily favor restraints (information-driven approach):
        // - VDW: 1.0 (baseline, prevents severe clashes)
        // - Elec: 0.2 (scaled for large values)
        // - Desolv: 1.0 (burial effects)
        // - AIR: 10.0 (DOMINANT - restraints are most important)
        let score = 1.0 * self.vdw + 0.2 * self.elec + 1.0 * self.desolv + 10.0 * self.air;

        self.fitness = score;
        self.fitness
    }
}
