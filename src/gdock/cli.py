# gdock
import argparse
import datetime
import logging
import sys

from gdock.modules.analysis import Analysis
from gdock.modules.ga import GeneticAlgorithm
from gdock.modules.geometry import Geometry
from gdock.modules.scoring import Scoring
from gdock.modules.initialize import Setup
from gdock.modules.structure import PDB, Restraint
from gdock.version import version
from gdock.modules.functions import random_quote

ga_log = logging.getLogger("ga_log")
ga_log.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter(
    " %(asctime)s %(module)s:%(lineno)d " "%(levelname)s - %(message)s"
)
ch.setFormatter(formatter)
ga_log.addHandler(ch)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("--dry", dest="dry", action="store_true", default=False)
    args = parser.parse_args()
    author, quote = random_quote()
    ga_log.info("==================================================================")
    ga_log.info("#                      Welcome to gdock!                          ")
    ga_log.info("#")
    ga_log.info(f'# "{quote}"')
    ga_log.info(f"#   -{author}")
    ga_log.info("#")
    ga_log.info("==================================================================")

    start_time = datetime.datetime.now()
    ga_log.info(f"Starting at {start_time.ctime()}")
    ga_log.info(f"Running {version}")

    ga_log.info("Setting up simulation")
    s = Setup(args.input_file)
    # TODO: implement setup file validation

    params = s.initialize()

    # 1. Load structure
    ga_log.info("Loading Structures module")
    input_molecules = PDB()
    input_molecules.load(params["mol_a"])
    input_molecules.load(params["mol_b"])

    # 2. Load restraints
    ga_log.info("Loading Restraints module")
    restraints = Restraint(input_molecules.raw_pdb)
    restraints.load(params["restraints_a"], "A")
    restraints.load(params["restraints_b"], "B")

    # 3. Position
    ga_log.info("Loading Geometry module")
    geo = Geometry(input_molecules, restraints)
    geo.calc_initial_position()
    initial_complex = geo.apply_transformation()

    # 4. Run GA
    ga_log.info("Loading Genetic Algorithm module")
    ga = GeneticAlgorithm(initial_complex, params)
    ga.setup()
    results = ga.run()

    if args.dry:
        ga_log.info("This is a dry run, stopping after models are generated")
        sys.exit()

    # 5. Scoring
    ga_log.info("Loading Scoring module")
    scoring = Scoring(results, params)
    results = scoring.score()

    # 6. Analysis
    ga_log.info("Loading Analysis module")
    ana = Analysis(results, params)
    ana.cluster()
    ana.evaluate()
    ana.output()

    s.clean()

    end_time = datetime.datetime.now()
    ga_log.info(f"Finishing at {end_time.ctime()}")
    duration = end_time - start_time
    ga_log.info(f"Took {duration}")

    ga_log.info("gdock complete.")
    # done (:


if __name__ == "__main__":
    main()
