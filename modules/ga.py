import os
import random
import numpy as np
import multiprocessing
from deap import base, creator, tools
from pyquaternion import Quaternion
from functions import calc_clash, timeit, dcomplex


class Population:
    def __init__(self, pioneer, target_chain, nproc):
        """

        :param pioneer:
        :param target_chain:
        :param nproc:
        """
        self.pioneer = pioneer
        self.nproc = nproc
        self.chain = target_chain
        if not os.path.isdir('pdbs/'):
            os.mkdir('pdbs/')

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
            for l in fh.readlines():
                if l.startswith('ATOM'):
                    chain = l[21]
                    x = float(l[31:38])
                    y = float(l[39:46])
                    z = float(l[47:54])
                    if chain not in pdb_dic:
                        pdb_dic[chain] = {'coord': [], 'raw': []}
                    pdb_dic[chain]['coord'].append((x,y,z))
                    pdb_dic[chain]['raw'].append(l)

        # rotate
        q = Quaternion(rotation)
        c = np.array(pdb_dic[target_chain]['coord'])
        center = c.mean(axis=0)
        c -= center
        r = np.array([q.rotate(e) for e in c])
        r += center
        pdb_dic[target_chain]['coord'] = list(r)

        # translate
        pass

        with open(output_fname, 'w') as out_fh:
            for chain in pdb_dic:
                for coord, line in zip(pdb_dic[chain]['coord'], pdb_dic[chain]['raw']):
                    new_x = f'{coord[0]:.3f}'.rjust(7, ' ')
                    new_y = f'{coord[1]:.3f}'.rjust(7, ' ')
                    new_z = f'{coord[2]:.3f}'.rjust(7, ' ')
                    new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}'
                    out_fh.write(new_line)
        out_fh.close()

        return output_fname

    @timeit
    def generate_pop(self, individuals):
        """

        :param individuals:
        """
        arg_list = []
        for i in individuals:
            output_fname = 'pdbs/gd_' + '_'.join(map(str, i)) + '.pdb'
            if not os.path.isfile(output_fname):
                arg_list.append((self.pioneer, i, self.chain, output_fname))
        # rotate(pdb_fname, rotation, target_chain, output_fname):
        pool = multiprocessing.Pool(processes=self.nproc)
        results = pool.starmap(self.rotate, arg_list)

    @staticmethod
    def output(coord_dic, output_fname):
        """

        :param coord_dic:
        :param output_fname:
        :return:
        """
        with open(output_fname, 'w') as out_fh:
            for chain in coord_dic:
                for coord, line in zip(coord_dic[chain]['coord'], coord_dic[chain]['raw']):
                    new_x = f'{coord[0]:.3f}'.rjust(7, ' ')
                    new_y = f'{coord[1]:.3f}'.rjust(7, ' ')
                    new_z = f'{coord[2]:.3f}'.rjust(7, ' ')
                    new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}'
                    out_fh.write(new_line)
        out_fh.close()
        return True

class GeneticAlgorithm(Population):

    def __init__(self, pioneer, population_size, number_of_generations, target_chain, nproc):
        super().__init__(pioneer, target_chain, nproc)
        # self.pop = Population(pioneer, target_chain)
        self.popsize = population_size
        self.ngen = number_of_generations
        self.generation_dic = {}

    @timeit
    def setup(self):
        """

        :return:
        """
        # -1 will optimize towards negative
        creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMax)

        # Individual and population functions
        toolbox = base.Toolbox()
        toolbox.register("attr_int", self.generate_individual, -1, 1)
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_int, n=4)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register("mutate", tools.mutPolynomialBounded, eta=0.8, low=-1, up=+1, indpb=0.2)
        toolbox.register("select", tools.selTournament, tournsize=2)
        toolbox.register("evaluate", self.fitness_function)

        pool = multiprocessing.Pool(processes=self.nproc)
        toolbox.register("map", pool.map)

        return toolbox

    def run(self, toolbox, cxpb=0.8, mutpb=0.2):
        """

        :param toolbox:
        :param cxpb:
        :param mutpb:
        :return:
        """
        pop = toolbox.population(n=self.popsize)
        for g in range(self.ngen):
            self.generation_dic[g] = {}
            offspring = toolbox.select(pop, len(pop))
            # Clone the population created
            offspring = list(map(toolbox.clone, offspring))

            # Apply crossover on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < cxpb:
                    toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            # Apply mutation on the offspring
            for mutant in offspring:
                if random.random() < mutpb:
                    toolbox.mutate(mutant)
                    del mutant.fitness.values

            # Generate the PDBs to be evaluated by the fitness function
            pop_pdbs = self.generate_pop(offspring)

            fitnesses = toolbox.map(toolbox.evaluate, offspring)
            for ind, fit in zip(offspring, fitnesses):
                ind.fitness.values = fit, fit

            # replace the old population by the offspring
            pop[:] = offspring

            for idx, ind in enumerate(pop):
                self.generation_dic[g][idx] = ind, []
                for fitness_v in ind.fitness.values:
                    self.generation_dic[g][idx][1].append(fitness_v)

            fitness_list = [self.generation_dic[g][f][1][0] for f in self.generation_dic[g]]
            print(f'+++ Generation {g}')
            print(f'Fitness: {fitness_list}')
            print(f'Average: {np.mean(fitness_list)}')
            print(f'Best: {min(fitness_list)}')


        return self.generation_dic

    @staticmethod
    def fitness_function(int_list):
        """

        :param int_list:
        :return:
        """
        rotated_pdb = 'pdbs/gd_' + '_'.join(map(str, int_list)) + '.pdb'
        # fit = calc_clash(rotated_pdb)
        fit = dcomplex(rotated_pdb)
        return fit

    @staticmethod
    def generate_individual(start, end):
        """

        :param start:
        :param end:
        :return:
        """
        # generate a random float between -1 and 1
        return round(random.choice(np.arange(start, end, 0.2)), 3)