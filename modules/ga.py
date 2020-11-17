import os
import secrets
import multiprocessing
import numpy as np
import toml as toml
from tempfile import NamedTemporaryFile
from deap import base, creator, tools
from scipy.spatial.transform import Rotation as R
from utils.files import get_full_path
from utils.functions import format_coords
from modules.fitness import run_foldx
import logging

ga_log = logging.getLogger('ga_log')

ga_params = toml.load(f"{get_full_path('etc')}/genetic_algorithm_params.toml")
secretsGenerator = secrets.SystemRandom()

# This needs to be outside..! https://github.com/rsteca/sklearn-deap/issues/59
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))  # -1 will optimize towards negative
creator.create("Individual", list, fitness=creator.FitnessMin)


class GeneticAlgorithm:

    def __init__(self, pioneer, run_params):
        """Initialize GeneticAlgorithm class."""
        self.run_params = run_params
        self.nproc = self.run_params['np']
        # self.ngen = ga_params['general']['number_of_generations']
        self.popsize = ga_params['general']['population_size']
        self.cxpb = ga_params['general']['crossover_probability']
        self.mutpb = ga_params['general']['mutation_probability']
        self.eta = ga_params['general']['eta']
        self.indpb = ga_params['general']['indpb']
        self.toolbox = None
        self.generation_dic = {}
        self.pioneer_dic = {}
        for line in pioneer.split('\n'):
            if line.startswith('ATOM'):
                chain = line[21]
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                if chain not in self.pioneer_dic:
                    self.pioneer_dic[chain] = {'coord': [], 'raw': []}
                self.pioneer_dic[chain]['coord'].append((x, y, z))
                self.pioneer_dic[chain]['raw'].append(line)

    # @timeit
    def setup(self):
        """Setup the genetic algorithm."""
        ga_log.debug('Creating the creator')

        # creator was here!

        # Individual and population functions
        ga_log.debug('Creating the individual and population functions')
        toolbox = base.Toolbox()
        toolbox.register("attr", self.generate_individual)
        toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.attr)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register("mutate_rot", tools.mutPolynomialBounded, eta=self.eta, low=0, up=360, indpb=self.indpb)
        toolbox.register("mutate_trans", tools.mutPolynomialBounded, eta=self.eta, low=-4, up=+4, indpb=self.indpb)
        toolbox.register("select", tools.selTournament, tournsize=2)
        toolbox.register("evaluate", self.fitness_function,
                         self.pioneer_dic)

        ga_log.debug('Creating the multiprocessing pool')
        pool = multiprocessing.Pool(processes=self.nproc)
        toolbox.register("map", pool.map)

        self.toolbox = toolbox

    def run(self):
        """Run the genetic algorithm."""
        ga_log.info('Running GA!')
        conv_l = []
        result_l = []
        run = True
        ga_log.info(f'Generations: Inf. Population: {self.popsize}')
        pop = self.toolbox.population(n=self.popsize)
        ngen = 1
        while run:
            self.generation_dic[ngen] = {}
            offspring = self.toolbox.select(pop, len(pop))
            # Clone the population created
            offspring = list(map(self.toolbox.clone, offspring))

            # Apply crossover on the offspring
            ga_log.debug('Applying crossover')
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if secretsGenerator.uniform(0, 1) < self.cxpb:
                    self.toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            # Apply mutation on the offspring
            ga_log.debug('Applying mutation')
            for mutant in offspring:
                if secretsGenerator.uniform(0, 1) < self.mutpb:
                    self.toolbox.mutate_rot(mutant[:3])
                    self.toolbox.mutate_trans(mutant[3:])
                    del mutant.fitness.values

            ga_log.debug('Calculating fitnessess')

            # =====
            #
            # uncoment below to debug the fitness function
            #
            # self.fitness_function(self.pioneer_dic, offspring[0])
            #
            # =====
            fitnesses = self.toolbox.map(self.toolbox.evaluate, offspring)
            for ind, fit in zip(offspring, fitnesses):
                ind.fitness.values = fit  # , fit

            # replace the old population by the offspring
            pop[:] = offspring

            for idx, ind in enumerate(pop):
                self.generation_dic[ngen][idx] = ind, []
                for fitness_v in ind.fitness.values:
                    self.generation_dic[ngen][idx][1].append(fitness_v)

            irmsd_list = [self.generation_dic[ngen][f][1][0] for f in self.generation_dic[ngen]]
            mean_fitness = np.mean(irmsd_list)
            std_fitness = np.std(irmsd_list)
            max_fitness = max(irmsd_list)
            min_fitness = min(irmsd_list)

            result_l.append(mean_fitness)
            if len(result_l) >= 2:
                conv = round(result_l[-1] - result_l[-2], 3)
            else:
                conv = round(.0, 3)

            conv_l.append(conv)

            ngen_str = str(ngen).rjust(3, '0')
            ga_log.info(f"Gen {ngen_str} fitness {mean_fitness:.2f} +- {std_fitness:.2f} [{max_fitness:.2f},"
                        f"{min_fitness:.2f}] ({conv:.3f})")

            if len(conv_l) >= 3 and sum(conv_l[-3:]) == .0:
                ga_log.info('Simulation converged, activating kill-switch!')
                run = False

            ngen += 1

        return self.generation_dic

    @staticmethod
    def fitness_function(pdb_dic, individual):
        """Calculate the fitness of an individual."""
        # use the chromossome and create the structure!
        c = np.array(pdb_dic['B']['coord'])
        # transform
        rot = R.from_euler('zyx', individual[:3])
        transl = individual[3:]

        center = c.mean(axis=0)
        c -= center
        c -= transl
        r = np.array([rot.apply(e) for e in c])
        r += center

        pdb_dic['B']['coord'] = list(r)

        # use a temporary file, nothing lasts forever
        pdb = NamedTemporaryFile(delete=False, suffix='.pdb')
        for chain in pdb_dic:
            for coord, line in zip(pdb_dic[chain]['coord'], pdb_dic[chain]['raw']):
                new_x, new_y, new_z = format_coords(coord)
                new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}\n'
                pdb.write(str.encode(new_line))
        pdb.close()

        # Calculate fitnesses!
        # ================================#
        energy = run_foldx(pdb.name)
        # ================================#

        # unlink the pdb so that it disappears
        os.unlink(pdb.name)

        # this must (?) be a list: github.com/DEAP/deap/issues/256
        return [energy]

    @staticmethod
    def generate_individual():
        """Generates the individual."""
        ind = [round(secretsGenerator.uniform(0, 360), 2),
               round(secretsGenerator.uniform(0, 360), 2),
               round(secretsGenerator.uniform(0, 360), 2),
               round(secretsGenerator.uniform(-4, 4), 3),
               round(secretsGenerator.uniform(-4, 4), 3),
               round(secretsGenerator.uniform(-4, 4), 3)]

        return ind
