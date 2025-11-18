use super::constants::MAX_DISPLACEMENT;
use crate::fitness;
use crate::restraints;
use crate::structure;
use crate::StdRng;
use core::f64::consts::PI;

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
        w_vdw: f64,
        w_elec: f64,
        w_desolv: f64,
        w_air: f64,
    ) -> f64 {
        let target_ligand = self.apply_genes(ligand);

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
        let score =
            w_vdw * self.vdw + w_elec * self.elec + w_desolv * self.desolv + w_air * self.air;

        self.fitness = score;
        self.fitness
    }
}
