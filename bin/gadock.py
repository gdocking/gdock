# GA needs to solve the Quaternion (grid position implemented later)
import argparse
import os
from modules.ga import GeneticAlgorithm
from modules.structure import PDB

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("input")
    parser.add_argument("npop", type=int)
    parser.add_argument("ngen", type=int)
    parser.add_argument("nproc", type=int)
    args = parser.parse_args()

    # 1. Load structure
    pdb = PDB()
    pdb.load(args.input)
    # 2. Position
    pass
    # DEV: randomize orientation of chain B
    pdb.randomize_rotation('B')
    pdb.output('rotated.pdb')
    # print(dcomplex('target-unbound.pdb'))

    # Proof-of-concept
    # from a unbound complex, randomize the ligand initial orientation and use GA to find it again, use clash as fitness
    ga = GeneticAlgorithm('rotated.pdb',
                          population_size=args.npop,
                          number_of_generations=args.ngen,
                          target_chain='B',
                          nproc = args.nproc)
    toolbox = ga.setup()
    result_dic = ga.run(toolbox)

    # make a nice movie!
    movie = []
    for gen in result_dic:
        best = []
        for ind in result_dic[gen]:
            seq, fitness = result_dic[gen][ind]
            name = 'pdbs/gd_' + '_'.join(map(str, seq)) + '.pdb'
            best.append((name, fitness[0]))
        best.sort(key=lambda x: x[1])
        middle_index = int(len(best) / 2)
        movie.append(best[middle_index][0])

    movie_str = ' '.join(movie)
    os.system(f'/Users/rodrigo/software/anaconda3/bin/pdb_mkensemble {movie_str} > movie.pdb')

    # done :)
