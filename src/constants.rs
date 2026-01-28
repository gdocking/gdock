// GA Parameters
pub const POPULATION_SIZE: u64 = 150;
pub const MUTATION_RATE: f64 = 0.1;
pub const CROSSOVER_RATE: f64 = 0.6;
pub const TOURNAMENT_SIZE: u64 = 3;
pub const ELITISM_COUNT: usize = 5; // Preserve top 5 individuals
pub const MAX_GENERATIONS: u64 = 250;
// NOTE: Changing the MAX_DISPLACEMENT from 15 to 50 greatly increased the score in debug from 0.7
// to 0.9!!
pub const MAX_DISPLACEMENT: f64 = 50.0;
pub const RANDOM_SEED: u64 = 42;

// Early stopping parameters
pub const ENABLE_EARLY_STOPPING: bool = true;
pub const CONVERGENCE_THRESHOLD: f64 = 0.001; // 0.1% improvement threshold
pub const CONVERGENCE_WINDOW: u64 = 10; // Stop if no improvement for this many generations

// Default weights
//  These defaults are based on a score calibration,
//   check the `calibration/` directory of the base
//   repository for more information
pub const DEFAULT_W_VDW: f64 = 0.4;
pub const DEFAULT_W_ELEC: f64 = 0.05;
pub const DEFAULT_W_DESOLV: f64 = 3.4;
pub const DEFAULT_W_AIR: f64 = 1.0;

/// Energy function weights for scoring
#[derive(Debug, Clone, Copy)]
pub struct EnergyWeights {
    pub vdw: f64,
    pub elec: f64,
    pub desolv: f64,
    pub air: f64,
}

impl Default for EnergyWeights {
    fn default() -> Self {
        Self {
            vdw: DEFAULT_W_VDW,
            elec: DEFAULT_W_ELEC,
            desolv: DEFAULT_W_DESOLV,
            air: DEFAULT_W_AIR,
        }
    }
}

impl EnergyWeights {
    pub fn new(vdw: f64, elec: f64, desolv: f64, air: f64) -> Self {
        Self {
            vdw,
            elec,
            desolv,
            air,
        }
    }
}

// Energy function cutoffs and parameters
pub const VDW_CUTOFF: f64 = 12.0; // Angstroms
pub const SOFTCORE_ALPHA: f64 = 0.5; // Softness parameter for soft-core LJ
pub const DESOLV_CUTOFF: f64 = 8.0; // Angstroms
pub const ELEC_CUTOFF: f64 = 15.0; // Angstroms
pub const ELEC_MIN_DISTANCE: f64 = 1.0; // Angstroms - avoid singularities
pub const AIR_FORCE_CONSTANT: f64 = 10.0; // kcal/mol per Å²

// VdW radii
// `https://en.wikipedia.org/wiki/Van_der_Waals_radius`
pub const HYDROGEN_RADIUS: f64 = 1.2;
pub const CARBON_RADIUS: f64 = 1.7;
pub const NITROGEN_RADIUS: f64 = 1.55;
pub const OXYGEN_RADIUS: f64 = 1.52;
pub const SULFUR_RADIUS: f64 = 1.8;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_energy_weights_new() {
        let weights = EnergyWeights::new(1.0, 0.5, 0.1, 100.0);
        assert_eq!(weights.vdw, 1.0);
        assert_eq!(weights.elec, 0.5);
        assert_eq!(weights.desolv, 0.1);
        assert_eq!(weights.air, 100.0);
    }

    #[test]
    fn test_energy_weights_default() {
        let weights = EnergyWeights::default();
        assert_eq!(weights.vdw, DEFAULT_W_VDW);
        assert_eq!(weights.elec, DEFAULT_W_ELEC);
        assert_eq!(weights.desolv, DEFAULT_W_DESOLV);
        assert_eq!(weights.air, DEFAULT_W_AIR);
    }

    #[test]
    fn test_energy_weights_custom_values() {
        let weights = EnergyWeights::new(2.0, 1.0, 0.2, 50.0);
        assert_eq!(weights.vdw, 2.0);
        assert_eq!(weights.elec, 1.0);
        assert_eq!(weights.desolv, 0.2);
        assert_eq!(weights.air, 50.0);
    }

    #[test]
    fn test_energy_weights_copy() {
        let weights1 = EnergyWeights::new(1.0, 0.5, 0.1, 100.0);
        let weights2 = weights1; // Copy
        assert_eq!(weights2.vdw, 1.0);
        assert_eq!(weights2.elec, 0.5);
    }

    #[test]
    fn test_energy_weights_clone() {
        let weights1 = EnergyWeights::new(1.0, 0.5, 0.1, 100.0);
        let weights2 = weights1; // EnergyWeights implements Copy
        assert_eq!(weights2.vdw, 1.0);
        assert_eq!(weights2.elec, 0.5);
    }
}
