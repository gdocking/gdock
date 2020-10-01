# GADock
import argparse

from modules.geometry import Geometry
from modules.setup import Setup
from modules.structure import PDB, Restraint
from modules.ga import GeneticAlgorithm
# from utils.functions import get_coords

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    args = parser.parse_args()

    s = Setup(args.input_file)
    if not s.validate():
        # raise GADockError
        exit()

    run_params = s.initialize()

    # 1. Load structure
    input_molecules = PDB()
    input_molecules.load(run_params['mol_a'])
    input_molecules.load(run_params['mol_b'])

    # 2. Load restraints
    restraints = Restraint(input_molecules.raw_pdb)
    restraints.load(run_params['restraints_a'], 'A')
    restraints.load(run_params['restraints_b'], 'B')

    # 3. Position
    geo = Geometry(input_molecules, restraints)
    geo.initial_position()

    # input_mol['A']
    # pass

    # DEV: randomize orientation of chain B
    # pdb.randomize_rotation('B')
    # pdb.output('init.pdb')

    # pdb.prepare_cell() # cell = grid

    # 3. Run GA
    ga = GeneticAlgorithm('examples/init.pdb',
                          target_chain='B',
                          nproc=args.np)
    toolbox = ga.setup()
    result_dic = ga.run(toolbox)
    output = ga.output('gadock.dat')
    # plot = ga.plot('plot.png')

    # done :)
