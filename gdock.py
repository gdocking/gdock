# gdock
import argparse
import datetime
import sys
from modules.setup import Setup
from modules.geometry import Geometry
from modules.structure import PDB, Restraint
from modules.ga import GeneticAlgorithm
from modules.scoring import Scoring
from modules.analysis import Analysis
from utils.functions import random_quote
from modules.version import CURRENT_VERSION
import logging

ga_log = logging.getLogger('ga_log')
ga_log.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter(' %(asctime)s %(module)s:%(lineno)d '
                              '%(levelname)s - %(message)s')
ch.setFormatter(formatter)
ga_log.addHandler(ch)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument('--dry', dest='dry', action='store_true', default=False)
    args = parser.parse_args()
    author, quote = random_quote()
    print('==================================================================')
    print('#                      Welcome to gdock!                          ')
    print('#')
    print(f'# "{quote}"')
    print(f'#   -{author}')
    print('#')
    print('==================================================================')

    start_time = datetime.datetime.now()
    ga_log.info(f'Starting at {start_time.ctime()}')
    ga_log.info(f'Running {CURRENT_VERSION}')

    ga_log.info('Setting up simulation')
    s = Setup(args.input_file)
    # TODO: implement setup file validation

    run_params, ga_params = s.initialize()

    # 1. Load structure
    ga_log.info('Loading Structures module')
    input_molecules = PDB()
    input_molecules.load(run_params['mol_a'])
    input_molecules.load(run_params['mol_b'])

    # 2. Load restraints
    ga_log.info('Loading Restraints module')
    restraints = Restraint(input_molecules.raw_pdb)
    restraints.load(run_params['restraints_a'], 'A')
    restraints.load(run_params['restraints_b'], 'B')

    # 3. Position
    ga_log.info('Loading Geometry module')
    geo = Geometry(input_molecules, restraints)
    geo.calc_initial_position()
    initial_complex = geo.apply_transformation()

    # 4. Run GA
    ga_log.info('Loading Genetic Algorithm module')
    ga = GeneticAlgorithm(initial_complex, run_params, ga_params)
    ga.setup()
    results = ga.run()

    if args.dry:
        ga_log.info('This is a dry run, stopping after models are generated')
        sys.exit()

    # 5. Scoring
    ga_log.info('Loading Scoring module')
    scoring = Scoring(results, run_params)
    results = scoring.score()

    # 6. Analysis
    ga_log.info('Loading Analysis module')
    ana = Analysis(results, run_params)
    ana.cluster(cutoff=0.6)
    ana.evaluate()
    ana.output()

    s.clean()

    end_time = datetime.datetime.now()
    ga_log.info(f'Finishing at {end_time.ctime()}')
    duration = end_time - start_time
    ga_log.info(f'Took {duration}')

    ga_log.info('gdock complete.')
    # done (:
