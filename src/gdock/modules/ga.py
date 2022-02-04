import copy
import logging
import multiprocessing
import os
import random
import sys
from tempfile import NamedTemporaryFile
from operator import attrgetter

import numpy as np
from deap import base, creator, tools

from gdock.modules.fitness import calc_satisfaction
from gdock.modules.geometry import rotate, translate
from gdock.modules.functions import format_coords, summary
from gdock.modules.error import GAError

ga_log = logging.getLogger("ga_log")


# This needs to be outside..! https://github.com/rsteca/sklearn-deap/issues/59
# -1 will optimize towards negative
# creator.create("FitnessMulti", base.Fitness, weights=(+1.0, -0.8))
creator.create("FitnessSingle", base.Fitness, weights=(+1.0,))
creator.create("Individual", list, fitness=creator.FitnessSingle)


# Custom GA functions
# selTournament


def select(individuals, k, tournsize):
    """Select the best individual"""
    # :param individuals: A list of individuals to select from.
    # :param k: The number of individuals to select.
    # :param tournsize: The number of individuals participating in each tournament.
    # :param fit_attr: The attribute of individuals to use as selection criterion
    # :returns: A list of selected individuals.

    # This is based on DEAP's selTournament but we account for
    #  the individuals being tuples
    chosen = []
    for i in range(k):
        # randomly pick 3
        aspirants = tools.selRandom(individuals, tournsize)
        # select the best of these 3
        sub_l = [(e, e[1].fitness) for e in aspirants]
        sub_l.sort(key=lambda x: x[1])
        best = sub_l[0][0]
        chosen.append(best)
    return chosen


def generate_individual():
    """Generates the individual."""
    ind = [
        random.randint(0, 360),  # nosec
        random.randint(0, 360),
        random.randint(0, 360),
        random.randint(-5, 5),
        random.randint(-5, 5),
        random.randint(-5, 5),
    ]

    return ind


def generate_population(func, structure_dictionary, pop_size):
    """Generate the population."""
    population = []

    total_complexes = len(structure_dictionary)
    if pop_size < total_complexes:
        raise GAError(
            f"Population size={pop_size} which is smaller then the number of pioneers={total_complexes}"
        )

    # each complex will be sampled a few times depending on the pop_size
    factor = int(pop_size / total_complexes)

    for pioneer in structure_dictionary:
        for _ in range(factor):
            population.append((pioneer, func()))

    return population


class GeneticAlgorithm:
    def __init__(self, pioneers, run_params, ga_params):
        """Initialize GeneticAlgorithm class."""
        self.random_seed = ga_params["parameters"]["random_seed"]
        self.run_params = run_params
        self.nproc = self.run_params["np"]
        self.structure_folder = f"{self.run_params['folder']}/structures"
        self.max_ngen = ga_params["general"]["max_number_of_generations"]
        self.popsize = ga_params["general"]["population_size"]
        self.cxpb = ga_params["general"]["crossover_probability"]
        self.mutpb = ga_params["general"]["mutation_probability"]
        self.eta = ga_params["general"]["eta"]
        self.indpb = ga_params["general"]["indpb"]
        self.conv_counter = ga_params["parameters"]["convergence_cutoff"]
        self.toolbox = None
        self.generation_dic = {}

        if self.popsize < len(pioneers):
            raise GAError(
                f"Population size={self.popsize} is smaller than the number of pioneers={len(pioneers)}"
            )
        factor = int(self.popsize / len(pioneers))
        self.pioneers = pioneers * factor
        # self.pioneer_dic = {}
        # for line in pioneer.split("\n"):
        #     if line.startswith("ATOM"):
        #         chain = line[21]
        #         x = float(line[31:38])
        #         y = float(line[39:46])
        #         z = float(line[47:54])
        #         if chain not in self.pioneer_dic:
        #             self.pioneer_dic[chain] = {"coord": [], "raw": []}
        #         self.pioneer_dic[chain]["coord"].append((x, y, z))
        #         self.pioneer_dic[chain]["raw"].append(line)
        # # TODO: find a better way of assigning chains here
        # self.pioneer_dic["A"]["restraints"] = run_params["restraints_a"]
        # self.pioneer_dic["B"]["restraints"] = run_params["restraints_b"]

    def setup(self):
        """Setup the genetic algorithm."""
        ga_log.debug("Creating the creator")

        # creator was here!

        # Individual and population functions
        ga_log.debug("Creating the individual and population functions")
        toolbox = base.Toolbox()
        toolbox.register("attr", generate_individual)
        toolbox.register(
            "individual", tools.initIterate, creator.Individual, toolbox.attr
        )
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register(
            "mutate_rot",
            tools.mutPolynomialBounded,
            eta=self.eta,
            low=0,
            up=360,
            indpb=self.indpb,
        )
        toolbox.register(
            "mutate_trans",
            tools.mutPolynomialBounded,
            eta=self.eta,
            low=-5,
            up=+5,
            indpb=self.indpb,
        )

        toolbox.register("select", tools.selTournament, tournsize=3)

        # toolbox.register("evaluate", self.fitness_function, self.pioneer_dic)
        toolbox.register("evaluate", self.fitness_function, self.pioneers)

        self.toolbox = toolbox

    def run(self):
        """Run the genetic algorithm."""
        ga_log.info("Running the Genetic Algorithm!")

        ga_log.debug("Creating the multiprocessing pool")
        pool = multiprocessing.Pool(processes=self.nproc)
        self.toolbox.register("map", pool.map)

        ga_log.info(f"Your random seed is: {self.random_seed}")
        random.seed(self.random_seed)
        variation_l = []
        result_l = []
        run = True
        ga_log.info(f"Population: {self.popsize} Max Generations:" f" {self.max_ngen}")
        ga_log.info("Gen ... mean +- sd (min, max)")
        pop = self.toolbox.population(n=self.popsize)
        ngen = 1
        while run:
            self.generation_dic[ngen] = {}

            offspring = self.toolbox.select(pop, len(pop))

            # Clone the population created
            offspring = list(map(self.toolbox.clone, offspring))

            # Apply crossover on the offspring
            ga_log.debug("Applying crossover")
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < self.cxpb:  # nosec
                    self.toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            # Apply mutation on the offspring
            ga_log.debug("Applying mutation")
            for mutant in offspring:
                if random.random() < self.mutpb:  # nosec
                    self.toolbox.mutate_rot(mutant[:3])
                    self.toolbox.mutate_trans(mutant[3:])
                    del mutant.fitness.values

            ga_log.debug("Calculating fitnessess")

            # =====
            #
            # uncoment below to debug the fitness function
            #
            # self.fitness_function(self.pioneers, (0, offspring[0]))
            #
            # =====
            # evaluate only the mutated and crossed
            # invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            invalid_ind = [
                (id, ind) for id, ind in enumerate(offspring) if not ind.fitness.valid
            ]
            try:
                fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind)
            except Exception as e:
                ga_log.warning(e)
                ga_log.error("Fitness function could not be executed")
                sys.exit()

            for ind, fit in zip(invalid_ind, fitnesses):
                ind[1].fitness.values = fit  # , fit

            # replace the old population by the offspring
            pop[:] = offspring

            satisfaction_values = []
            for idx, ind in enumerate(pop):
                fitness_v = ind.fitness.values
                self.generation_dic[ngen][idx] = {
                    "individual": ind,
                    "fitness": fitness_v,
                }
                # energy = ind.fitness.values[1]
                # energy_values.append(energy)

                satisfaction = ind.fitness.values[0]
                satisfaction_values.append(satisfaction)

            # energy_summary = summary(energy_values)
            satisfaction_summary = summary(satisfaction_values)

            # result_l.append(energy_summary['mean'])
            result_l.append(satisfaction_summary["mean"])

            last_results = result_l[-self.conv_counter :]
            variation = min(last_results) - max(last_results)

            variation_l.append(variation)

            ngen_str = str(ngen).rjust(3, "0")
            # ga_log.info(f"Gen {ngen_str} ({satisfaction_summary['mean']:.2f}"
            #             f", {energy_summary['mean']:.3f})")
            ga_log.info(
                f"Gen {ngen_str} {satisfaction_summary['mean']:.2f} "
                f"+- {satisfaction_summary['std']:.2f} "
                f"({satisfaction_summary['min']:.2f}, "
                f"{satisfaction_summary['max']:.2f}) "
                f"eval={len(invalid_ind)}"
            )

            if ngen == self.max_ngen:
                ga_log.info(
                    f"Simulation reached maximum number of "
                    f"generations, stopping at {ngen}."
                )
                run = False

            if len(result_l) >= self.conv_counter:
                convergence = []
                for var in variation_l[-self.conv_counter :]:
                    if abs(var) < 0.01:
                        convergence.append(True)
                    else:
                        convergence.append(False)

                if all(convergence):
                    ga_log.info('Simulation "converged"')
                    ga_log.info(
                        f"Absolute mean fitness variation is < .01 for"
                        f" last {self.conv_counter} generations"
                    )
                    ga_log.info(f"Stopped at generation {ngen}")
                    run = False

            ngen += 1

        pool.close()
        pool.join()

        # ga complete, generate the epoch
        self._generate_epoch()

        return self.generation_dic

    @staticmethod
    def fitness_function(pioneers, element):
        """Calculate the fitness of an individual."""
        # fun fact: if you do not use deepcopy here, the fitness will depend on the number of processors
        #  since we are using a shared data structure
        ident, individual = element
        _pioneers = copy.deepcopy(pioneers)
        pioneer = _pioneers[ident - 1]
        individual_dic = {}
        for line in pioneer.structure.split("\n"):
            if line.startswith("ATOM"):
                chain = line[21]
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                if chain not in individual_dic:
                    individual_dic[chain] = {"coord": [], "raw": []}
                individual_dic[chain]["coord"].append((x, y, z))
                individual_dic[chain]["raw"].append(line)

        individual_str = " ".join([f"{j:.2f}" for j in individual])
        c = np.array(individual_dic["B"]["coord"])

        translation_center = individual[3:]
        rotation_angles = individual[:3]

        translated_coords = translate(c, translation_center)
        rotated_coords = rotate(translated_coords, rotation_angles)

        individual_dic["B"]["coord"] = list(rotated_coords)

        # use a temporary file to keep the execution simple with
        #  some I/O trade-off
        pdb = NamedTemporaryFile(delete=False, suffix=".pdb")
        for chain in individual_dic:
            coord_l = individual_dic[chain]["coord"]
            raw_l = individual_dic[chain]["raw"]
            for coord, line in zip(coord_l, raw_l):
                new_x, new_y, new_z = format_coords(coord)
                new_line = (
                    f"{line[:30]} {new_x} {new_y} {new_z} " f"{line[55:]}" + os.linesep
                )
                pdb.write(str.encode(new_line))
        pdb.close()

        # Calculate fitnesses!
        # ================================#
        # energy = run_dcomplex(pdb.name)
        satisfaction = calc_satisfaction(
            pdb.name,
            pioneer.restraint_i_resnums,
            pioneer.restraint_j_resnums,
            # individual_dic["A"]["restraints"],
            # individual_dic["B"]["restraints"],
        )
        # ================================#

        # unlink the pdb so that it disappears
        os.unlink(pdb.name)

        # return [satisfaction, energy]
        ga_log.debug(f"{individual_str} {satisfaction:.2f} {pdb.name}")
        # this must (?) be a list: github.com/DEAP/deap/issues/256
        return [satisfaction]

    def _generate_epoch(self):
        """Read information from simulation and generate structures."""
        ga_log.info("Dumping structures for this epoch")

        pool = multiprocessing.Pool(processes=self.nproc)

        done_dic = {}
        clone_counter = 0
        for gen in self.generation_dic:
            for idx in self.generation_dic[gen]:
                individual = self.generation_dic[gen][idx]["individual"]
                individual_as_str = "_".join(map(str, individual))
                if individual_as_str not in done_dic:
                    pdb_name = (
                        f"{self.structure_folder}/"
                        f"{str(gen).rjust(4, '0')}"
                        f"_{str(idx).rjust(4, '0')}.pdb"
                    )
                    done_dic[individual_as_str] = pdb_name
                    self.generation_dic[gen][idx]["structure"] = pdb_name
                    self.generation_dic[gen][idx]["clone"] = None

                    pool.apply_async(
                        self._recreate, args=(self.pioneers[idx], individual, pdb_name)
                    )
                else:
                    # this is a clone, discard it
                    clone_name = done_dic[individual_as_str]
                    self.generation_dic[gen][idx]["clone"] = clone_name
                    self.generation_dic[gen][idx]["structure"] = None
                    clone_counter += 1

        pool.close()
        pool.join()

        ga_log.debug(f"{clone_counter} clones were found")

    @staticmethod
    def _recreate(complex, individual, pdb_name):
        """Use the chromossome information and create the structure."""
        # input_dic = copy.deepcopy(complex)
        _complex = copy.deepcopy(complex)
        # c = np.array(input_dic["B"]["coord"])
        c = _complex.coords_j

        translation_center = individual[3:]
        rotation_angles = individual[:3]

        translated_coords = translate(c, translation_center)
        rotated_coords = rotate(translated_coords, rotation_angles)

        # input_dic["B"]["coord"] = list(rotated_coords)
        _complex.coords_j = list(rotated_coords)

        with open(pdb_name, "w") as fh:
            for coord, line in zip(_complex.coords_i, _complex.pdb_i):
                new_x, new_y, new_z = format_coords(coord)
                new_line = (
                    f"{line[:30]} {new_x} {new_y} {new_z}" f" {line[55:]}" + os.linesep
                )
                fh.write(new_line)
            for coord, line in zip(_complex.coords_j, _complex.pdb_j):
                new_x, new_y, new_z = format_coords(coord)
                new_line = (
                    f"{line[:30]} {new_x} {new_y} {new_z}" f" {line[55:]}" + os.linesep
                )
                fh.write(new_line)
