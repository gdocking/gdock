import copy
import logging
import multiprocessing
import os
import random
import sys
from tempfile import NamedTemporaryFile

import numpy as np
from deap import base, creator, tools

from gdock.modules.fitness import calc_satisfaction, calc_haddock_score, run_dcomplex
from gdock.modules.geometry import Geometry
from gdock.modules.functions import format_coords, summary, tidy

ga_log = logging.getLogger("ga_log")


# This needs to be outside..! https://github.com/rsteca/sklearn-deap/issues/59
# optimize restraints to max and energy to min
creator.create("FitnessMulti", base.Fitness, weights=(1.0, -0.8))
creator.create("Individual", list, fitness=creator.FitnessMulti)

IND_DB = {}
ROT_RANGE = 360
TRANS_RANGE = 2


class GeneticAlgorithm:
    def __init__(self, pioneer, params):
        """Initialize GeneticAlgorithm class."""
        self.params = params
        self.structure_folder = f"{params['folder']}/structures"
        self.random_seed = params["main"]["random_seed"]
        self.nproc = params["main"]["number_of_processors"]
        self.energy_function = params["main"]["scoring_function"]

        self.max_ngen = params["ga"]["max_number_of_generations"]
        self.popsize = params["ga"]["population_size"]
        self.cxpb = params["ga"]["crossover_probability"]
        self.mutpb = params["ga"]["mutation_probability"]
        self.eta = params["ga"]["eta"]
        self.indpb = params["ga"]["indpb"]
        self.conv_counter = params["ga"]["convergence_counter"]

        self.toolbox = None
        self.generation_dic = {}
        self.diff_dic = {"satisfaction": [], "energy": []}

        self.pioneer_dic = self._load_pioneer(pioneer)

    def _load_pioneer(self, pioneer):
        """Load the Pioneer into the data structure."""
        pioneer_dic = {}
        for line in pioneer.split("\n"):
            if line.startswith("ATOM"):
                chain = line[21]
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                if chain not in pioneer_dic:
                    pioneer_dic[chain] = {"coord": [], "raw": []}
                pioneer_dic[chain]["coord"].append((x, y, z))
                pioneer_dic[chain]["raw"].append(line)
        # TODO: find a better way of assigning chains here
        pioneer_dic["A"]["restraints"] = self.params["restraints_a"]
        pioneer_dic["B"]["restraints"] = self.params["restraints_b"]

        return pioneer_dic

    def setup(self):
        """Setup the genetic algorithm."""
        ga_log.debug("Creating the creator")

        # creator was here!

        # Individual and population functions
        ga_log.debug("Creating the individual and population functions")
        toolbox = base.Toolbox()
        toolbox.register("attr", self.generate_individual)
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
            up=ROT_RANGE,
            indpb=self.indpb,
        )
        toolbox.register(
            "mutate_trans",
            tools.mutPolynomialBounded,
            eta=self.eta,
            low=-TRANS_RANGE,
            up=+TRANS_RANGE,
            indpb=self.indpb,
        )

        toolbox.register("select", tools.selTournament, tournsize=3)

        toolbox.register(
            "evaluate", self.fitness_function, self.pioneer_dic, self.energy_function
        )

        self.toolbox = toolbox

    def run(self):
        """Run the genetic algorithm."""
        ga_log.debug("Creating the multiprocessing pool")
        pool = multiprocessing.Pool(processes=self.nproc)
        self.toolbox.register("map", pool.map)

        ga_log.info("#" * 72)

        ga_log.info("Genetic Algorithm parameters: ")
        ga_log.info(f"  + population_size: {self.popsize}")
        ga_log.info(f"  + max_number_of_generations: {self.max_ngen}")
        ga_log.info(f"  + crossover_probability: {self.cxpb}")
        ga_log.info(f"  + mutation_probability: {self.mutpb}")
        ga_log.info(f"  + eta (crowding degree of mut): {self.eta}")
        ga_log.info(f"  + indpb (independent prob): {self.indpb}")
        ga_log.info(f"  + convergence_counter: {self.conv_counter}")
        ga_log.info("General parameters: ")
        ga_log.info(f"  + random_seed: {self.random_seed}")
        ga_log.info(f"  + scoring_function: {self.energy_function}")
        ga_log.info(f"  + nproc (number of processors): {self.nproc}")

        random.seed(self.random_seed)

        ga_log.info("#" * 72)

        ga_log.info("> satisfaction = ratio of satisfied restraints")

        if self.energy_function == "dcomplex":
            ga_log.info("> energy = dcomplex_energy * -1")

        elif self.energy_function == "haddock":
            ga_log.info("> energy = haddock_score")

        ga_log.info("#" * 72)
        ga_log.info("Gen ... |  satisfaction  |      energy     | variation")

        run = True
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
            # uncomment below to debug the fitness function
            #
            # self.fitness_function(self.pioneer_dic, offspring[0])
            #
            # =====
            # evaluate only the mutated and crossed
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            try:
                fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind)
            except Exception as e:
                ga_log.warning(e)
                ga_log.error("Fitness function could not be executed")
                sys.exit()

            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit  # , fit

            # replace the old population by the offspring
            pop[:] = offspring

            for idx, ind in enumerate(pop):
                satisfaction, energy = ind.fitness.values
                self.generation_dic[ngen][idx] = {
                    "individual": ind,
                    "satisfaction": satisfaction,
                    "energy": energy,
                }

            satisfaction_diff, energy_diff = self.check_diff(generation=ngen)

            satisfaction_summary = summary(
                [
                    self.generation_dic[ngen][_ind]["satisfaction"]
                    for _ind in self.generation_dic[ngen]
                ]
            )
            energy_summary = summary(
                [
                    self.generation_dic[ngen][_ind]["energy"]
                    for _ind in self.generation_dic[ngen]
                ]
            )
            ngen_str = str(ngen).rjust(3, "0")

            ga_log.info(
                f"Gen {ngen_str} "
                f"| {satisfaction_summary['mean']:.3f} "
                f"+- {satisfaction_summary['std']:.3f} "
                f"| {energy_summary['mean']:.2f} "
                f"+- {energy_summary['std']:.2f} "
                f"| {satisfaction_diff:.3f}, {energy_diff:.3f} "
            )

            if self.convergence():
                ga_log.info("#" * 72)
                ga_log.info(
                    f'Simulation "converged", no change over the latest {self.conv_counter} generations'
                )
                ga_log.info(f"Stopped at generation {ngen}")
                run = False

            if ngen == self.max_ngen:
                ga_log.info("#" * 72)
                ga_log.info(
                    f"Simulation reached maximum number of "
                    f"generations, stopping at {ngen}."
                )
                run = False

            ngen += 1

        pool.close()
        pool.join()

        # ga complete, generate the epoch
        self._generate_epoch()

        return self.generation_dic

    @staticmethod
    def fitness_function(pdb_dic, energy_function, individual):
        """Calculate the fitness of an individual."""
        # use the chromossome and create the structure!

        # fun fact: if you do not use deepcopy here,
        #  the fitness will depend on the number of processors
        #  since pdb_dic is a shared data structure
        individual_str = " ".join([f"{j:.2f}" for j in individual])
        if individual_str in IND_DB:
            return IND_DB[individual_str]

        individual_dic = copy.deepcopy(pdb_dic)
        c = np.array(individual_dic["B"]["coord"])

        rotation_angles = individual[:3]
        translation_center = individual[3:]

        translated_coords = Geometry.translate(c, translation_center)
        rotated_coords = Geometry.rotate(translated_coords, rotation_angles)

        individual_dic["B"]["coord"] = list(rotated_coords)

        # use a temporary file to keep the execution simple with
        #  some I/O trade-off
        pdb = NamedTemporaryFile(delete=False, suffix=".pdb", mode="w")
        pdb_str = ""
        for chain in individual_dic:
            coord_l = individual_dic[chain]["coord"]
            raw_l = individual_dic[chain]["raw"]
            for coord, line in zip(coord_l, raw_l):
                new_x, new_y, new_z = format_coords(coord)
                new_line = (
                    f"{line[:30]} {new_x} {new_y} {new_z} " f"{line[55:]}" + os.linesep
                )
                pdb_str += new_line

        # Tidy and write the PDB
        pdb.write(tidy(pdb_str))

        # Calculate fitnesses!
        # ================================#
        satisfaction = calc_satisfaction(
            pdb.name,
            individual_dic["A"]["restraints"],
            individual_dic["B"]["restraints"],
        )

        if energy_function == "dcomplex":
            dcomplex_score = run_dcomplex(pdb.name)
            # in dcomplex the higher energy the better, -1 here to flip the direction
            #  and optimize towards the lowest
            energy_score = dcomplex_score * -1

        elif energy_function == "haddock":
            energy_score = calc_haddock_score(pdb.name)

        # ================================#

        # unlink the pdb so that it disappears
        os.unlink(pdb.name)

        ga_log.debug(
            f"{individual_str} {satisfaction:.2f} {energy_score:.2f} {pdb.name}"
        )

        IND_DB[individual_str] = satisfaction, energy_score

        # this must (?) be a list: github.com/DEAP/deap/issues/256
        return [satisfaction, energy_score]

    @staticmethod
    def generate_individual():
        """Generates the individual."""
        ind = [
            random.uniform(0, ROT_RANGE),  # nosec
            random.uniform(0, ROT_RANGE),
            random.uniform(0, ROT_RANGE),
            random.uniform(-TRANS_RANGE, +TRANS_RANGE),
            random.uniform(-TRANS_RANGE, +TRANS_RANGE),
            random.uniform(-TRANS_RANGE, +TRANS_RANGE),
        ]

        return [round(n, 3) for n in ind]

    def convergence(self):
        """Check if simulation has converged."""
        # Conversion means no change in the last N generations
        latest_satisfaction_l = self.diff_dic["satisfaction"][-self.conv_counter :]
        latest_energy_l = self.diff_dic["energy"][-self.conv_counter :]

        assert len(latest_satisfaction_l) == len(latest_energy_l)

        if (
            len(latest_satisfaction_l) == self.conv_counter
            and all(v == 0.0 for v in latest_satisfaction_l)
            and all(v == 0.0 for v in latest_energy_l)
        ):
            return True
        else:
            return False

    def check_diff(self, generation):
        """Check if the values increased or not."""

        if generation == 1:
            avg_energy_diff = 0.0
            avg_satisfaction_diff = 0.0
        else:
            current_gendic = self.generation_dic[generation]
            previous_gendic = self.generation_dic[generation - 1]

            current_avg_satisfaction = np.average(
                [current_gendic[ind]["satisfaction"] for ind in current_gendic]
            )
            previous_avg_satisfaction = np.average(
                [previous_gendic[ind]["satisfaction"] for ind in previous_gendic]
            )
            avg_satisfaction_diff = current_avg_satisfaction - previous_avg_satisfaction

            current_avg_energy = np.average(
                [current_gendic[ind]["energy"] for ind in current_gendic]
            )
            previous_avg_energy = np.average(
                [previous_gendic[ind]["energy"] for ind in previous_gendic]
            )
            avg_energy_diff = current_avg_energy - previous_avg_energy

        self.diff_dic["satisfaction"].append(avg_satisfaction_diff)
        self.diff_dic["energy"].append(avg_energy_diff)

        return avg_satisfaction_diff, avg_energy_diff

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
                        self._recreate, args=(self.pioneer_dic, individual, pdb_name)
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
    def _recreate(input_structure_dic, individual, pdb_name):
        """Use the chromossome information and create the structure."""
        input_dic = copy.deepcopy(input_structure_dic)
        c = np.array(input_dic["B"]["coord"])

        translation_center = individual[3:]
        rotation_angles = individual[:3]

        translated_coords = Geometry.translate(c, translation_center)
        rotated_coords = Geometry.rotate(translated_coords, rotation_angles)

        input_dic["B"]["coord"] = list(rotated_coords)

        with open(pdb_name, "w") as fh:
            for chain in input_dic:
                coord_l = input_dic[chain]["coord"]
                raw_l = input_dic[chain]["raw"]
                for coord, line in zip(coord_l, raw_l):
                    new_x, new_y, new_z = format_coords(coord)
                    new_line = (
                        f"{line[:30]} {new_x} {new_y} {new_z}"
                        f" {line[55:]}" + os.linesep
                    )
                    fh.write(new_line)
        fh.close()
