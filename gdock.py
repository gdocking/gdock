# gdock
import argparse
from modules.setup import Setup
from modules.geometry import Geometry
from modules.structure import PDB, Restraint
from modules.ga import GeneticAlgorithm
from modules.analysis import Analysis
from utils.functions import random_quote
from modules.version import CURRENT_VERSION
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
    author, quote = random_quote()
    print('============================================================================')
    print('#                            Welcome to gdock!                              ')
    print('#')
    print(f'# "{quote}"')
    print(f'#   -{author}')
    print('#')
    print('============================================================================')
    ga_log.info(f'Running {CURRENT_VERSION}')

    ga_log.info('Setting up simulation')
    s = Setup(args.input_file)
    # TODO: implement setup file validation

    run_params, ga_params = s.initialize()

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

    ga = GeneticAlgorithm(initial_complex, run_params, ga_params)
    ga.setup()
    results = ga.run()

    # 5. Analysis
    ga_log.info('Loading Analysis')
    ana = Analysis(initial_complex, results, run_params)
    ana.generate_structures()
    ana.cluster()
    ana.evaluate()
    ana.output()

    s.clean()

    ga_log.info('gdock finished.')
    # done (:
