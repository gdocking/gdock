# GADock
import argparse
import os
import shutil
from modules.ga import GeneticAlgorithm
from modules.structure import PDB

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("input")
    parser.add_argument("--np", type=int, default=4)
    args = parser.parse_args()

    # 1. Load structure
    pdb = PDB()
    pdb.load(args.input)
    # 2. Position
    pass
    # DEV: randomize orientation of chain B
    pdb.randomize_rotation('B')
    pdb.output('rotated.pdb')

    # Proof-of-concept
    # from a unbound complex, randomize the ligand initial orientation and use GA to find it again, use clash as fitness
    ga = GeneticAlgorithm('rotated.pdb',
                          target_chain='B',
                          nproc=args.np)
    toolbox = ga.setup()
    result_dic = ga.run(toolbox)
    output = ga.output('gadock.dat')
    plot = ga.plot('plot.png')

    # done :)
