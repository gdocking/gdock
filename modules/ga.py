import os
import secrets
import multiprocessing
import numpy as np
import toml as toml
from tempfile import NamedTemporaryFile
from deap import base, creator, tools
from pyquaternion import Quaternion
from utils.files import get_full_path
from utils.functions import format_coords
from modules.fitness import calc_irmsd
import logging

ga_log = logging.getLogger('ga_log')

ga_params = toml.load(f"{get_full_path('etc')}/genetic_algorithm_params.toml")
secretsGenerator = secrets.SystemRandom()


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

        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))  # -1 will optimize towards negative
        creator.create("Individual", list, fitness=creator.FitnessMin)

        # Individual and population functions
        ga_log.debug('Creating the individual and population functions')
        toolbox = base.Toolbox()
        toolbox.register("attr", self.generate_individual)
        toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.attr)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register("mutate_rot", tools.mutPolynomialBounded, eta=self.eta, low=-1, up=+1, indpb=self.indpb)
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
        kill_counter = 0
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
                    self.toolbox.mutate_rot(mutant[:4])
                    self.toolbox.mutate_trans(mutant[4:])
                    del mutant.fitness.values

            ga_log.debug('Calculating fitnessess')
            # self.fitness_function(self.pioneer_dic, offspring[0])
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

            ga_log.info(f" Gen {ngen}: irmsd: {mean_fitness:.2f} ± {std_fitness:.2f} [{max_fitness:.2f},"
                        f" {min_fitness:.2f}] ({conv:.3f})")

            if len(conv_l) >= 5 and sum(conv_l[-5:]) == .0:
                ga_log.info('Simulation converged, activating kill-switch!')
                run = False

            ngen += 1

        return self.generation_dic

    def output(self):
        """Output the fitness and individual properties."""
        folder = self.run_params['folder']
        output_f = f'{folder}/gdock.dat'
        ga_log.info(f'Writing output to {output_f}')
        with open(output_f, 'w') as fh:
            for gen in self.generation_dic:
                for ind in self.generation_dic[gen]:
                    individual = ' '.join(map(str, self.generation_dic[gen][ind][0]))
                    fitness = self.generation_dic[gen][ind][1][0]
                    tbw = f'{gen},{ind},{fitness},{individual}\n'
                    fh.write(tbw)
        fh.close()

    # def plot(self, plot_name):
    #     """
    #
    #     :type plot_name: string
    #     """
    #     import pandas as pd
    #     import seaborn as sns
    #     import matplotlib.pyplot as plt
    #     sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    #
    #     def label(color, plot_label):
    #         ax = plt.gca()
    #         ax.text(0, .2, plot_label, fontweight="bold", color=color,
    #                 ha="left", va="center", transform=ax.transAxes)
    #
    #     result_l = []
    #     for gen in self.generation_dic:
    #         for ind in self.generation_dic[gen]:
    #             name = 'pdbs/gd_' + '_'.join(map("{:.2f}".format, self.generation_dic[gen][ind][0])) + '.pdb'
    #             fitness = self.generation_dic[gen][ind][1][0]
    #             result_l.append((gen, ind, fitness, name))
    #     df = pd.DataFrame(result_l)
    #     df.columns = ['generation', 'individual', 'fitness', 'pdb']
    #
    #     pal = sns.cubehelix_palette(n_colors=len(df.generation.unique()), rot=0, light=.7, reverse=False)
    #     g = sns.FacetGrid(df, row="generation", hue="generation", aspect=15, height=.5, palette=pal)
    #     g.map(sns.kdeplot, "fitness", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
    #     g.map(sns.kdeplot, "fitness", clip_on=False, color="w", lw=2, bw=.2)
    #     g.map(plt.axhline, y=0, lw=2, clip_on=False)
    #
    #     g.map(label, "fitness")
    #     # Set the subplots to overlap
    #     g.fig.subplots_adjust(hspace=-.25)
    #
    #     # Remove axes details that don't play well with overlap
    #     g.set_titles("")
    #     g.set(yticks=[])
    #     g.despine(bottom=True, left=True)
    #     plt.savefig(plot_name)

    @staticmethod
    def fitness_function(pdb_dic, individual):
        """Calculate the fitness of an individual."""
        # use the chromossome and create the structure!
        c = np.array(pdb_dic['B']['coord'])
        # transform
        rot_q = Quaternion(individual[:4])
        transl = individual[4:]

        center = c.mean(axis=0)
        c -= center
        c -= transl
        r = np.array([rot_q.rotate(e) for e in c])
        r += center

        pdb_dic['B']['coord'] = list(r)

        # use a temporary file, nothing lasts forever
        pdb = NamedTemporaryFile(delete=False)
        for chain in pdb_dic:
            for coord, line in zip(pdb_dic[chain]['coord'], pdb_dic[chain]['raw']):
                new_x, new_y, new_z = format_coords(coord)
                new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}\n'
                pdb.write(str.encode(new_line))
        pdb.close()

        # Calculate fitnesses!
        # ================================#
        irmsd = calc_irmsd(pdb.name)
        # clash = calc_clash(rotated_pdb)
        # center_distance = calc_centerdistance(rotated_pdb)
        # fit = dcomplex(rotated_pdb)
        # ================================#

        # unlink the pdb so that it disappears
        os.unlink(pdb.name)

        # this must (?) be a list: github.com/DEAP/deap/issues/256
        return [irmsd]

    @staticmethod
    def generate_individual():
        """Generates the individual."""
        ind = [round(secretsGenerator.uniform(-1, 1), 2),
               round(secretsGenerator.uniform(-1, 1), 2),
               round(secretsGenerator.uniform(-1, 1), 2),
               round(secretsGenerator.uniform(-1, 1), 2),
               round(secretsGenerator.uniform(-4, 4), 3),
               round(secretsGenerator.uniform(-4, 4), 3),
               round(secretsGenerator.uniform(-4, 4), 3)]

        return ind
