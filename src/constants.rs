pub const POPULATION_SIZE: u64 = 250;
pub const MUTATION_RATE: f64 = 0.1;
pub const CROSSOVER_RATE: f64 = 0.6;
pub const TOURNAMENT_SIZE: u64 = 3;
pub const ELITISM_COUNT: usize = 5; // Preserve top 5 individuals (2%)
pub const MAX_GENERATIONS: u64 = 250;
pub const MAX_DISPLACEMENT: f64 = 15.0; // Keep molecules within interaction range!
pub const RANDOM_SEED: u64 = 42;

// Early stopping parameters
pub const ENABLE_EARLY_STOPPING: bool = true;
pub const CONVERGENCE_THRESHOLD: f64 = 0.001; // 0.1% improvement threshold
pub const CONVERGENCE_WINDOW: u64 = 20; // Stop if no improvement for this many generations

// VdW radii
// `https://en.wikipedia.org/wiki/Van_der_Waals_radius`
pub const HYDROGEN_RADIUS: f64 = 1.2;
pub const CARBON_RADIUS: f64 = 1.7;
pub const NITROGEN_RADIUS: f64 = 1.55;
pub const OXYGEN_RADIUS: f64 = 1.52;
pub const SULFUR_RADIUS: f64 = 1.8;
