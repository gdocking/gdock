# GADock
import argparse
from modules.structure import PDB
from modules.ga import GeneticAlgorithm

if __name__ == '__main__':

    # Proof-of-concept
    # from a unbound complex, randomize the ligand initial orientation and use GA to find it again

    parser = argparse.ArgumentParser(description='')
    # parser.add_argument("input")
    parser.add_argument("--np", type=int, default=4)
    args = parser.parse_args()

    # 1. Load structure
    # pdb = PDB()
    # pdb.load(args.input)

    # 2. Position
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
