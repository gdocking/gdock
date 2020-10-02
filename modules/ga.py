import glob
import os
import random
import numpy as np
import multiprocessing
import toml as toml
from deap import base, creator, tools
from dual_quaternions import DualQuaternion
from pyquaternion import Quaternion
from utils.files import get_full_path
from utils.functions import format_coords
from modules.fitness import calc_irmsd, calc_clash, calc_centerdistance
import logging
ga_log = logging.getLogger('ga_log')

ga_params = toml.load(f"{get_full_path('etc')}/genetic_algorithm_params.toml")


class Population:
    def __init__(self, pioneer, run_params):
        """

        :param pioneer:
        :param run_params:
        """
        self.run_params = run_params
        self.pioneer = f"{self.run_params['folder']}/begin/pioneer.pdb"
        with open(self.pioneer, 'w') as fh:
            ga_log.debug(f'Saving pioneer to {self.pioneer}')
            fh.write(pioneer)
        self.nproc = self.run_params['np']
        self.chain = 'B'

    @staticmethod
    def rotate(pdb_fname, rotation, target_chain, output_fname):
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

        # TODO: This here is the space exploration function
        #  Implement something for the translations!

        # # transform
        # dq = DualQuaternion.from_dq_array(q)
        # rot_q = Quaternion(dq.quat_pose_array()[:4])
        # transl = dq.quat_pose_array()[4:]
        # c = np.array(pdb_dic[target_chain]['coord'])
        # center = c.mean(axis=0)
        # c -= center
        # r = np.array([rot_q.rotate(e) for e in c])
        # r -= transl
        # r += center
        # pdb_dic[target_chain]['coord'] = list(r)

        # rotate
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
    def generate_pop(self, individuals, clean=True):
        """

        :param individuals:
        :param clean:
        """
        if clean:
            # Delete individuals from previous generations
            files = glob.glob('gen/*pdb')
            for f in files:
                os.remove(f)

        arg_list = []
        for i in individuals:
            output_fname = self.run_params['folder'] + '/gen/gd_' + '_'.join(map("{:.2f}".format, i)) + '.pdb'
            if not os.path.isfile(output_fname):
                arg_list.append((self.pioneer, i, self.chain, output_fname))
        # self.rotate(self.pioneer, i, self.chain, output_fname)
        pool = multiprocessing.Pool(processes=self.nproc)
        pool.starmap(self.rotate, arg_list)

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

    def __init__(self, pioneer, run_params):
        """

        :param pioneer:
        :param target_chain:
        :param nproc:
        """
        super().__init__(pioneer, run_params)
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

        creator.create("FitnessMax", base.Fitness, weights=(-1.0,)) #  -1 will optimize towards negative
        creator.create("Individual", list, fitness=creator.FitnessMax)

        # keep this for the future, to apply multiple fitnesses =====
        # creator.create("FitnessMulti", base.Fitness, weights=(-0.8, -0.5, -1.0))
        # creator.create("Individual", list, fitness=creator.FitnessMulti)
        # =========

        # Individual and population functions
        ga_log.debug('Creating the individual and population functions')
        toolbox = base.Toolbox()
        toolbox.register("attr_int", self.generate_individual, -1, 1)
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_int, n=4)
        # toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_int, n=8)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register("mutate", tools.mutPolynomialBounded, eta=self.eta, low=-1, up=+1, indpb=self.indpb)
        toolbox.register("select", tools.selTournament, tournsize=2)
        toolbox.register("evaluate", self.fitness_function)

        ga_log.debug('Creating the multiprocessing pool')
        pool = multiprocessing.Pool(processes=self.nproc)
        toolbox.register("map", pool.map)

        return toolbox

    def run(self, toolbox):
        """

        :param toolbox: DEAP toolbox object
        :return: dictionary {generation: [individual, (score1, score2, ...)]}
        """
        pop = toolbox.population(n=self.popsize)
        for g in range(self.ngen):
            self.generation_dic[g] = {}
            offspring = toolbox.select(pop, len(pop))
            # Clone the population created
            offspring = list(map(toolbox.clone, offspring))

            # Apply crossover on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < self.cxpb:
                    toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            # Apply mutation on the offspring
            for mutant in offspring:
                if random.random() < self.mutpb:
                    toolbox.mutate(mutant)
                    del mutant.fitness.values

            # Generate the PDBs to be evaluated by the fitness function
            self.generate_pop(offspring, clean=False)

            fitnesses = toolbox.map(toolbox.evaluate, offspring)
            for ind, fit in zip(offspring, fitnesses):
                ind.fitness.values = fit  # , fit

            # replace the old population by the offspring
            pop[:] = offspring

            for idx, ind in enumerate(pop):
                self.generation_dic[g][idx] = ind, []
                for fitness_v in ind.fitness.values:
                    self.generation_dic[g][idx][1].append(fitness_v)

            irmsd_list = [self.generation_dic[g][f][1][0] for f in self.generation_dic[g]]
            clash_list = [self.generation_dic[g][f][1][1] for f in self.generation_dic[g]]
            dista_list = [self.generation_dic[g][f][1][2] for f in self.generation_dic[g]]

            print("+" * 42)
            print(f"+ Generation {str(g).rjust(2, '0')} ++++++++++++++++++++++++++")
            print("+" * 42)
            print(f"  irmsd: {np.mean(irmsd_list):.2f} ± {np.std(irmsd_list):.2f} [{max(irmsd_list):.2f},"
                  f" {min(irmsd_list):.2f}]")
            print(f"  clash: {np.mean(clash_list):.2f} ± {np.std(clash_list):.2f} [{max(clash_list):.2f},"
                  f" {min(clash_list):.2f}]")
            print(f"  cdist: {np.mean(dista_list):.2f} ± {np.std(dista_list):.2f} [{max(dista_list):.2f},"
                  f" {min(dista_list):.2f}]")

        # save only the last generation for the future
        self.generate_pop(pop, clean=False)

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
        # return irmsd, clash, center_distance

        return irmsd

    @staticmethod
    def generate_individual(start, end):
        """

        :param start:
        :param end:
        :return:
        """
        # generate a random float between -1 and 1
        return round(random.choice(np.arange(start, end, 0.1)), 3)
