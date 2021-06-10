import random
import argparse
import logging
import multiprocessing
from deap import base, creator, tools
import scipy.stats
import pandas as pd


opt_log = logging.getLogger('opt_log')
opt_log.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter(' %(asctime)s L%(lineno)d '
                              '%(levelname)s > %(message)s')
ch.setFormatter(formatter)
opt_log.addHandler(ch)

# Create the creator #
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)


class OptimizationGA:

    def __init__(self, data_f, max_gen, pop, nproc=1, cxpb=0.5, mutpb=0.2):
        self.data = pd.read_table(data_f)
        self.max_gen = max_gen
        self.pop = pop
        self.cxpb = cxpb
        self.mutpb = mutpb
        self.nproc = nproc

    def setup(self):
        """Setup the optimization GA."""
        opt_log.info('Setting up GA')
        toolbox = base.Toolbox()
        toolbox.register('attr', self.generate_individual)

        # Register the toolbox #
        toolbox.register("individual",
                         tools.initIterate,
                         creator.Individual,
                         toolbox.attr)

        toolbox.register("population",
                         tools.initRepeat,
                         list,
                         toolbox.individual)

        toolbox.register("evaluate",
                         self.calc_correlation,
                         self.data)

        toolbox.register("mate", tools.cxTwoPoint)

        toolbox.register("mutate",
                         tools.mutPolynomialBounded,
                         eta=0.1,
                         low=-1,
                         up=+1,
                         indpb=0.4)

        toolbox.register("select", tools.selTournament, tournsize=3)

        self.toolbox = toolbox

    def run(self):
        """Run the optimization."""
        pool = multiprocessing.Pool(processes=self.nproc)
        self.toolbox.register("map", pool.map)

        random.seed(42)
        pop = self.toolbox.population(n=self.pop)

        opt_log.info('Starting optimization')

        # fitnesses = self.toolbox.map(self.toolbox.evaluate, pop)
        # for ind, fit in zip(pop, fitnesses):
        #     ind.fitness.values = fit

        # fits = [ind.fitness.values[0] for ind in pop]
        fits = [.0]

        # Begin the evolution
        opt_log.info('Starting the evolution')
        opt_log.info(f'Max generations={self.max_gen} Pop={len(pop)}')
        gen = 0
        while max(fits) < 1 and gen < self.max_gen:
            # A new generation
            gen += 1

            # Select the next generation individuals
            offspring = self.toolbox.select(pop, len(pop))

            # Clone the selected individuals
            offspring = list(map(self.toolbox.clone, offspring))

            # Apply crossover and mutation on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):

                # cross two individuals with probability CXPB
                if random.random() < self.cxpb:
                    self.toolbox.mate(child1, child2)

                    # fitness values of the children
                    # must be recalculated later
                    del child1.fitness.values
                    del child2.fitness.values

            for mutant in offspring:

                # mutate an individual with probability MUTPB
                if random.random() < self.mutpb:
                    self.toolbox.mutate(mutant)
                    # https://stackoverflow.com/a/44722548
                    del mutant.fitness.values

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            opt_log.debug(f"Evaluated {len(invalid_ind)} individuals")

            # The population is entirely replaced by the offspring
            pop[:] = offspring

            # Gather all the fitnesses in one list and print the stats
            fits = [ind.fitness.values[0] for ind in pop]

            length = len(pop)
            mean = sum(fits) / length
            sum2 = sum(x * x for x in fits)
            std = abs(sum2 / length - mean**2)**0.5

            opt_log.info(f"Gen {gen} {mean:.2f} +- {std:.2f} ({min(fits):.2f},"
                         f"{max(fits):.2f})")

        pool.close()
        pool.join()

        opt_log.info("Optimization complete")

        best_ind = tools.selBest(pop, 1)[0]
        opt_log.info(f"Optimal scoring function is: (energy * "
                     f"{best_ind[0]:.2f}) / (satisfaction *"
                     f" {best_ind[1]:.2f})")
        opt_log.info(f"This function will give a correlation of"
                     f" {best_ind.fitness.values[0]:.2f}")

    @staticmethod
    def generate_individual():
        """Generate random weights."""
        return [random.uniform(-1, +1),
                random.uniform(-1, +1)]

    @staticmethod
    def calc_correlation(data, individual):
        """Calculate the correlation between score and irmsd."""
        w_0, w_1 = individual
        score_l = []
        for e, s in zip(data['energy'], data['fitness']):
            gdock_score = (float(e) * w_0) / (float(s) * w_1)
            score_l.append(gdock_score)

        x = data['irmsd']
        y = pd.array(score_l)
        corr, p = scipy.stats.pearsonr(x, y)

        return corr,


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("data_f")
    args = parser.parse_args()
    # data_f = ('/Users/rodrigo/repos/gdock/benchmark/'
    #           'benchmark_v1.1.0_nocluster.dat')
    ga = OptimizationGA(args.data_f, max_gen=1000, pop=300, nproc=4)
    ga.setup()
    ga.run()
