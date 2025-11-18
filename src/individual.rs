#[derive(Debug, Clone)]
pub struct Population {
    pub chromosomes: Vec<Chromosome>,
    pub receptor: Molecule,
    pub ligand: Molecule,
    pub reference: Molecule,
    pub generation: u64,
}

impl Population {
    pub fn new(
        individuals: Vec<Chromosome>,
        receptor: Molecule,
        ligand: Molecule,
        reference: Molecule,
    ) -> Population {
        Population {
            chromosomes: individuals,
            receptor,
            ligand,
            reference,
            generation: 0,
        }
    }

    pub fn eval_fitness(&mut self) {
        // Calculate the fitness in parallel
        self.chromosomes.par_iter_mut().for_each(|chromossome| {
            // TODO: Move this whole fitness evaluation to a method inside Chromosome
            let mol = chromossome.apply_genes(self.ligand.clone());

            // Calculate the Fitness

            // let clashes = fitness::calc_clashes(&mol, &self.receptor);
            // let vdw_e = fitness::vdw_energy(&mol, &self.receptor);
            // let air_e = fitness::calculate_eair(&mol, &self.receptor);

            // chromossome.fitness = vdw_e + air_e;

            // ------
            // Calculate the metrics
            let complex = metrics::Complex::new(&self.receptor, &mol, &self.reference);
            // let (rmsd, fnat) = complex.;

            // ------
            chromossome.rmsd = complex.rmsd();
            chromossome.fnat = complex.fnat();
            chromossome.fitness = chromossome.rmsd;

            // println!("Evaluating fitness for chromossome {:?}", chromossome);

            // chromossome.model = mol;

            // // Write the fitness to a file
            // let mut output_file =
            //     std::fs::File::create(format!("dbg/{}_{}.fit", chromossome.id, self.generation))
            //         .unwrap();
            // write!(output_file, "{}", chromossome.fitness).unwrap();

            // // Write the PDB
            // let output_fname = format!("dbg/{}_{}.pdb", chromossome.id, self.generation);
            // write_pdb(&chromossome.model, &output_fname);
            // panic!();
        });
    }

    /// Evolve the population
    pub fn evolve(&self, rng: &mut StdRng) -> Population {
        // Tournament selection
        let mut new_population = self.tournament_selection(rng);

        // Crossover population
        // TODO: Improve this crossover part, we could do a weighted crossover based on the fitness
        //  Currently we are just doing a window of 2 elements
        new_population.chromosomes = new_population
            .chromosomes
            .windows(2) // Creates overlapping slices of 2 elements
            .map(|window| crossover(rng, &window[0], &window[1]))
            .collect();

        // Mutatation
        new_population
            .chromosomes
            .iter_mut()
            .for_each(|individual| {
                individual.mutate(rng, constants::MUTATION_RATE);
            });

        new_population
    }

    /// Tournament selection
    fn tournament_selection(&self, rng: &mut StdRng) -> Population {
        let mut offspring = Population::new(
            Vec::with_capacity(self.size()),
            self.receptor.clone(),
            self.ligand.clone(),
            self.reference.clone(),
        );

        offspring.generation = self.generation + 1;

        for _ in 0..self.size() {
            let tournament_individuals = (0..TOURNAMENT_SIZE)
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
