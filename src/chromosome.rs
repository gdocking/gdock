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
    // -------------------------------
    // This is just for debugging, remove it later
    // pub model: Molecule,
    pub rmsd: f64,
    pub fnat: f64,
    // -------------------------------
}

impl Chromosome {
    pub fn new(rng: &mut StdRng) -> Chromosome {
        // let mut rng = rand::thread_rng();

        let alpha = rng.gen_range(0.0..=2.0 * PI);
        let beta = rng.gen_range(0.0..=2.0 * PI);
        let gamma = rng.gen_range(0.0..=2.0 * PI);

        let displacement_x = rng.gen_range(-MAX_DISPLACEMENT..=MAX_DISPLACEMENT);
        let displacement_y = rng.gen_range(-MAX_DISPLACEMENT..=MAX_DISPLACEMENT);
        let displacement_z = rng.gen_range(-MAX_DISPLACEMENT..=MAX_DISPLACEMENT);

        // let genes =
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
            // model: Molecule::new(),
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

        let vdw_e = fitness::vdw_energy(receptor, &target_ligand);
        let elec_e = fitness::coulombic_energy(receptor, &target_ligand);
        let restraints_ratio = fitness::satisfaction_ratio(restraints, receptor, &target_ligand);

        let score = -1.0 * elec_e + 0.2 * vdw_e - 0.8 * restraints_ratio;

        self.fitness = score;
        self.fitness
    }
}
