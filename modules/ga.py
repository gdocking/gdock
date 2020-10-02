import glob
import os
import random
import numpy as np
import multiprocessing
import toml as toml
from deap import base, creator, tools
from pyquaternion import Quaternion
from utils.files import get_full_path
from utils.functions import format_coords
from modules.fitness import calc_irmsd, calc_clash, calc_centerdistance
import logging
ga_log = logging.getLogger('ga_log')

ga_params = toml.load(f"{get_full_path('etc')}/genetic_algorithm_params.toml")


class Population:
    def __init__(self, pioneer, grid, run_params):
        """

        :param pioneer:
        :param run_params:
        """
        self.run_params = run_params
        self.grid = grid
        self.pioneer = f"{self.run_params['folder']}/begin/pioneer.pdb"
        with open(self.pioneer, 'w') as fh:
            ga_log.debug(f'Saving pioneer to {self.pioneer}')
            fh.write(pioneer)
        self.nproc = self.run_params['np']
        self.chain = 'B'

    @staticmethod
    def explore(pdb_fname, individual, target_chain, output_fname, grid):
        """

        :param pdb_fname:
        :param rotation:
        :param target_chain:
        :param output_fname:
        :return:
        """
        pdb_dic = {}
        with open(pdb_fname, 'r') as fh:
            for line in fh.readlines():
                if line.startswith('ATOM'):
                    chain = line[21]
                    x = float(line[31:38])
                    y = float(line[39:46])
                    z = float(line[47:54])
                    if chain not in pdb_dic:
                        pdb_dic[chain] = {'coord': [], 'raw': []}
                    pdb_dic[chain]['coord'].append((x, y, z))
                    pdb_dic[chain]['raw'].append(line)

        # translate
        # TEMPORARY: Randomly move it to another grid point
        #  this makes no sense and is here just for implementation purposes
        c = np.array(pdb_dic[target_chain]['coord'])
        trans_chance = individual[-1]
        if trans_chance >= .1:
            grid_point, grid_coord = random.choice(list(grid.items()))
            c -= grid_coord

        # rotate
        rotation = individual[:4]
        q = Quaternion(rotation)
        c = np.array(pdb_dic[target_chain]['coord'])
        center = c.mean(axis=0)
        c -= center
        r = np.array([q.rotate(e) for e in c])
        r += center
        pdb_dic[target_chain]['coord'] = list(r)

        with open(output_fname, 'w') as out_fh:
            for chain in pdb_dic:
                for coord, line in zip(pdb_dic[chain]['coord'], pdb_dic[chain]['raw']):
                    new_x, new_y, new_z = format_coords(coord)
                    new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}'
                    out_fh.write(new_line)
        out_fh.close()

        return output_fname

    # @timeit
    def generate_pop(self, individuals):
        """

        :param individuals:
        :param clean:
        """

        arg_list = []
        for i in individuals:
            output_fname = self.run_params['folder'] + '/gen/gd_' + '_'.join(map("{:.2f}".format, i)) + '.pdb'
            if not os.path.isfile(output_fname):
                arg_list.append((self.pioneer, i, self.chain, output_fname, self.grid))
        # self.explore(self.pioneer, i, self.chain, output_fname, self.grid)
        ga_log.debug(f'Population generation len(arg_list): {len(arg_list)} nproc: {self.nproc}')
        pool = multiprocessing.Pool(processes=self.nproc)
        pool.starmap(self.explore, arg_list)

    # @staticmethod
    # def output(coord_dic, output_fname):
    #     """
    #
    #     :param coord_dic:
    #     :param output_fname:
    #     :return:
    #     """
    #     with open(output_fname, 'w') as out_fh:
    #         for chain in coord_dic:
    #             for coord, line in zip(coord_dic[chain]['coord'], coord_dic[chain]['raw']):
    #                 new_x, new_y, new_z = format_coords(coord)
    #                 new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}'
    #                 out_fh.write(new_line)
    #     out_fh.close()
    #     return True


class GeneticAlgorithm(Population):

    def __init__(self, pioneer, ligand_grid, run_params):
        """

        :param pioneer:
        :param target_chain:
        :param nproc:
        """
        super().__init__(pioneer, ligand_grid, run_params)
        # self.pop = Population(pioneer, target_chain)
        self.ngen = ga_params['general']['number_of_generations']
        self.popsize = ga_params['general']['population_size']
        self.cxpb = ga_params['general']['crossover_probability']
        self.mutpb = ga_params['general']['mutation_probability']
        self.eta = ga_params['general']['eta']
        self.indpb = ga_params['general']['indpb']
        self.generation_dic = {}

    # @timeit
    def setup(self):
        """

        :return:
        """
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
        toolbox.register("mutate", tools.mutPolynomialBounded, eta=self.eta, low=-1, up=+1, indpb=self.indpb)
        toolbox.register("select", tools.selTournament, tournsize=2)
        toolbox.register("evaluate", self.fitness_function)

        ga_log.debug('Creating the multiprocessing pool')
        pool = multiprocessing.Pool(processes=self.nproc)
        toolbox.register("map", pool.map)

        self.toolbox = toolbox

    def run(self):
        """

        :param toolbox: DEAP toolbox object
        :return: dictionary {generation: [individual, (score1, score2, ...)]}
        """
        ga_log.info('Running GA!')
        ga_log.info(f'Generations: {self.ngen} Population: {self.popsize}')
        pop = self.toolbox.population(n=self.popsize)
        for g in range(self.ngen):
            self.generation_dic[g] = {}
            offspring = self.toolbox.select(pop, len(pop))
            # Clone the population created
            offspring = list(map(self.toolbox.clone, offspring))

            # Apply crossover on the offspring
            ga_log.debug('Applying crossover')
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < self.cxpb:
                    self.toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            # Apply mutation on the offspring
            ga_log.debug('Applying mutation')
            for mutant in offspring:
                if random.random() < self.mutpb:
                    self.toolbox.mutate(mutant)
                    del mutant.fitness.values

            # Generate the PDBs to be evaluated by the fitness function
            ga_log.debug('Generating population')
            self.generate_pop(offspring)

            ga_log.debug('Calculating fitnessess')
            fitnesses = self.toolbox.map(self.toolbox.evaluate, offspring)
            for ind, fit in zip(offspring, fitnesses):
                ind.fitness.values = fit  # , fit

            # replace the old population by the offspring
            pop[:] = offspring

            for idx, ind in enumerate(pop):
                self.generation_dic[g][idx] = ind, []
                for fitness_v in ind.fitness.values:
                    self.generation_dic[g][idx][1].append(fitness_v)

            irmsd_list = [self.generation_dic[g][f][1][0] for f in self.generation_dic[g]]

            ga_log.info(f" Gen {g}: irmsd: {np.mean(irmsd_list):.2f} ± {np.std(irmsd_list):.2f} [{max(irmsd_list):.2f},"
                  f" {min(irmsd_list):.2f}]")

        # save only the last generation for the future
        self.generate_pop(pop)

        return self.generation_dic

    def output(self, output_f):
        with open(output_f, 'w') as fh:
            for gen in self.generation_dic:
                for ind in self.generation_dic[gen]:
                    name = 'pdbs/gd_' + '_'.join(map("{:.2f}".format, self.generation_dic[gen][ind][0])) + '.pdb'
                    fitness = self.generation_dic[gen][ind][1][0]
                    tbw = f'{gen},{ind},{fitness},{name}\n'
                    fh.write(tbw)
        fh.close()
        return True

    def plot(self, plot_name):
        """

        :type plot_name: string
        """
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

        def label(color, plot_label):
            ax = plt.gca()
            ax.text(0, .2, plot_label, fontweight="bold", color=color,
                    ha="left", va="center", transform=ax.transAxes)

        result_l = []
        for gen in self.generation_dic:
            for ind in self.generation_dic[gen]:
                name = 'pdbs/gd_' + '_'.join(map("{:.2f}".format, self.generation_dic[gen][ind][0])) + '.pdb'
                fitness = self.generation_dic[gen][ind][1][0]
                result_l.append((gen, ind, fitness, name))
        df = pd.DataFrame(result_l)
        df.columns = ['generation', 'individual', 'fitness', 'pdb']

        pal = sns.cubehelix_palette(n_colors=len(df.generation.unique()), rot=0, light=.7, reverse=False)
        g = sns.FacetGrid(df, row="generation", hue="generation", aspect=15, height=.5, palette=pal)
        g.map(sns.kdeplot, "fitness", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
        g.map(sns.kdeplot, "fitness", clip_on=False, color="w", lw=2, bw=.2)
        g.map(plt.axhline, y=0, lw=2, clip_on=False)

        g.map(label, "fitness")
        # Set the subplots to overlap
        g.fig.subplots_adjust(hspace=-.25)

        # Remove axes details that don't play well with overlap
        g.set_titles("")
        g.set(yticks=[])
        g.despine(bottom=True, left=True)
        plt.savefig(plot_name)

    @staticmethod
    def fitness_function(int_list):
        """

        :param int_list:
        :return:
        """
        # FIXME: Need to also give the PATH to this function, somehow
        rotated_pdb = '/Users/rodrigo/repos/gadock/dev/gen/gd_' + '_'.join(map("{:.2f}".format, int_list)) + '.pdb'
        irmsd = calc_irmsd(rotated_pdb)
        # clash = calc_clash(rotated_pdb)
        # center_distance = calc_centerdistance(rotated_pdb)
        # fit = dcomplex(rotated_pdb)

        # this MUST (?) be a list
        # github.com/DEAP/deap/issues/256
        return [irmsd]

    @staticmethod
    def generate_individual():
        """

        :param start:
        :param end:
        :return:
        """

        ind = []
        ind.append(round(random.choice(np.arange(-1, +1, 0.1)), 3))
        ind.append(round(random.choice(np.arange(-1, +1, 0.1)), 3))
        ind.append(round(random.choice(np.arange(-1, +1, 0.1)), 3))
        ind.append(round(random.choice(np.arange(-1, +1, 0.1)), 3))
        ind.append(random.random())

        return ind
