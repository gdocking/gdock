//! Tracks diverse high-quality solutions across evolution generations.
use std::f64::consts::PI;

use crate::chromosome::Chromosome;
use crate::constants::{
    HALL_OF_FAME_MAX_SIZE, HALL_OF_FAME_TOP_K, HOF_UNIQUENESS_ROTATION_THRESHOLD,
    HOF_UNIQUENESS_TRANSLATION_THRESHOLD,
};

#[derive(Debug, Clone)]
pub struct HallOfFameEntry {
    pub genes: [f64; 6],
    pub fitness: f64,
    pub vdw: f64,
    pub elec: f64,
    pub desolv: f64,
    pub air: f64,
    pub clash: f64,
}

#[derive(Debug)]
pub struct HallOfFame {
    entries: Vec<HallOfFameEntry>,
    max_size: usize,
}

impl HallOfFame {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            max_size: HALL_OF_FAME_MAX_SIZE,
        }
    }

    pub fn with_capacity(max_size: usize) -> Self {
        Self {
            entries: Vec::new(),
            max_size,
        }
    }

    /// Try to add an entry if it is unique enough. Returns `true` if added.
    #[allow(clippy::too_many_arguments)]
    pub fn try_add(
        &mut self,
        genes: &[f64],
        fitness: f64,
        vdw: f64,
        elec: f64,
        desolv: f64,
        air: f64,
        clash: f64,
    ) -> bool {
        if genes.len() != 6 {
            return false;
        }
        let new_genes: [f64; 6] = [genes[0], genes[1], genes[2], genes[3], genes[4], genes[5]];
        if !self.is_unique(&new_genes) {
            return false;
        }
        self.entries.push(HallOfFameEntry {
            genes: new_genes,
            fitness,
            vdw,
            elec,
            desolv,
            air,
            clash,
        });
        if self.entries.len() > self.max_size {
            self.prune();
        }
        true
    }

    fn is_unique(&self, new_genes: &[f64; 6]) -> bool {
        !self
            .entries
            .iter()
            .any(|e| Self::genes_are_similar(new_genes, &e.genes))
    }

    fn genes_are_similar(a: &[f64; 6], b: &[f64; 6]) -> bool {
        for i in 0..3 {
            if Self::angular_difference(a[i], b[i]) > HOF_UNIQUENESS_ROTATION_THRESHOLD {
                return false;
            }
        }
        for i in 3..6 {
            if (a[i] - b[i]).abs() > HOF_UNIQUENESS_TRANSLATION_THRESHOLD {
                return false;
            }
        }
        true
    }

    fn angular_difference(a: f64, b: f64) -> f64 {
        let diff = (a - b).abs();
        if diff > PI {
            2.0 * PI - diff
        } else {
            diff
        }
    }

    fn prune(&mut self) {
        self.entries
            .sort_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap());
        self.entries.truncate(self.max_size);
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn entries(&self) -> &[HallOfFameEntry] {
        &self.entries
    }

    /// Add the top-K unique chromosomes from a population snapshot.
    pub fn add_from_population(&mut self, chromosomes: &[Chromosome]) {
        let mut indices: Vec<usize> = (0..chromosomes.len()).collect();
        indices.sort_by(|&a, &b| {
            chromosomes[a]
                .fitness
                .partial_cmp(&chromosomes[b].fitness)
                .unwrap()
        });
        let mut added = 0;
        for &idx in &indices {
            if added >= HALL_OF_FAME_TOP_K {
                break;
            }
            let chr = &chromosomes[idx];
            if self.try_add(
                &chr.genes,
                chr.fitness,
                chr.vdw,
                chr.elec,
                chr.desolv,
                chr.air,
                chr.clash,
            ) {
                added += 1;
            }
        }
    }
}

impl Default for HallOfFame {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_hall_of_fame_new() {
        let hof = HallOfFame::new();
        assert!(hof.is_empty());
        assert_eq!(hof.len(), 0);
    }

    #[test]
    fn test_hall_of_fame_add_entry() {
        let mut hof = HallOfFame::new();
        let genes = [0.0_f64; 6];
        assert!(hof.try_add(&genes, -100.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        assert_eq!(hof.len(), 1);
    }

    #[test]
    fn test_hall_of_fame_rejects_duplicates() {
        let mut hof = HallOfFame::new();
        let genes = [0.0_f64; 6];
        hof.try_add(&genes, -100.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert!(!hof.try_add(&genes, -100.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        assert_eq!(hof.len(), 1);
    }

    #[test]
    fn test_hall_of_fame_accepts_different_genes() {
        let mut hof = HallOfFame::new();
        let genes1 = [0.0_f64; 6];
        let genes2 = [PI, PI, PI, 0.0, 0.0, 0.0];
        hof.try_add(&genes1, -100.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert!(hof.try_add(&genes2, -90.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        assert_eq!(hof.len(), 2);
    }

    #[test]
    fn test_hall_of_fame_accepts_different_translations() {
        let mut hof = HallOfFame::new();
        let genes1 = [0.0_f64; 6];
        let genes2 = [0.0, 0.0, 0.0, 10.0, 10.0, 10.0];
        hof.try_add(&genes1, -100.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert!(hof.try_add(&genes2, -90.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        assert_eq!(hof.len(), 2);
    }

    #[test]
    fn test_hall_of_fame_rejects_similar_genes() {
        let mut hof = HallOfFame::new();
        let genes1 = [0.0_f64; 6];
        let genes2 = [0.1, 0.1, 0.1, 0.5, 0.5, 0.5];
        hof.try_add(&genes1, -100.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert!(!hof.try_add(&genes2, -90.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        assert_eq!(hof.len(), 1);
    }

    #[test]
    fn test_hall_of_fame_prune_keeps_best() {
        let mut hof = HallOfFame::with_capacity(3);
        for i in 0..5 {
            let genes = [i as f64, 0.0, 0.0, (i as f64) * 10.0, 0.0, 0.0];
            hof.try_add(&genes, -(100.0 - i as f64 * 10.0), 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        assert_eq!(hof.len(), 3);
        let entries = hof.entries();
        assert!(entries[0].fitness <= entries[1].fitness);
        assert!(entries[1].fitness <= entries[2].fitness);
    }

    #[test]
    fn test_angular_difference_same() {
        assert!((HallOfFame::angular_difference(0.0, 0.0) - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_angular_difference_wrap_around() {
        let diff = HallOfFame::angular_difference(0.0, 2.0 * PI - 0.1);
        assert!(diff < 0.2, "Should handle wrap-around");
    }

    #[test]
    fn test_angular_difference_opposite() {
        let diff = HallOfFame::angular_difference(0.0, PI);
        assert!((diff - PI).abs() < 0.001);
    }
}
