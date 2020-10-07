# GADock
import argparse
from modules.geometry import Geometry
from modules.setup import Setup
from modules.structure import PDB, Restraint
from modules.ga import GeneticAlgorithm
import logging

ga_log = logging.getLogger('ga_log')
ga_log.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter(' %(asctime)s %(module)s:%(lineno)d %(levelname)s - %(message)s')
ch.setFormatter(formatter)
ga_log.addHandler(ch)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    args = parser.parse_args()

    ga_log.info('Setting up simulation')
    s = Setup(args.input_file)
    # TODO: implement setup file validation

    run_params = s.initialize()

    # 1. Load structure
    ga_log.info('Loading structures')
    input_molecules = PDB()
    input_molecules.load(run_params['mol_a'])
    input_molecules.load(run_params['mol_b'])

    # 2. Load restraints
    ga_log.info('Loading restraints')
    restraints = Restraint(input_molecules.raw_pdb)
    restraints.load(run_params['restraints_a'], 'A')
    restraints.load(run_params['restraints_b'], 'B')

    # 3. Position
    ga_log.info('Loading geometry')
    geo = Geometry(input_molecules, restraints)
    geo.calc_initial_position()
    initial_complex = geo.apply_transformation()

    # 4. Run GA
    ga_log.info('Loading Genetic Algorithm')
    ga = GeneticAlgorithm(initial_complex, run_params)
    ga.setup()
    result_dic = ga.run()
    output = ga.output()
    # plot = ga.plot('plot.png')

    # done :)
