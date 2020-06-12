# read the output file and create a trajectory based on the fitness
import os
from operator import itemgetter
from utils.functions import add_dummy_center

movie_dic = {}
with open('gadock.dat','r') as fh:
    for l in fh.readlines():
        data = l.split(',')
        gen = data[0]
        ind = data[1]
        fitness = float(data[2])
        pdb = data[3].split('\n')[0]
        if gen not in movie_dic:
            movie_dic[gen] = []
        movie_dic[gen].append((pdb, fitness))

movie_l = []
for g in movie_dic:
    # select worst, middle, best
    gen = sorted(movie_dic[g], key=itemgetter(1))
    best = gen[0][0]
    middle = gen[int((len(gen)- 1)/2)][0]
    worst = gen[-1][0]
    movie_l.append(best)
    movie_l.append(middle)
    movie_l.append(worst)

for i, pdb in enumerate(movie_l):
    if not os.path.isfile(pdb):
        print(f'File not found {pdb}')
    else:
        c_pdb = add_dummy_center(pdb)
        os.system(f"/Users/rodrigo/software/anaconda3/bin/pdb_tidy {c_pdb} > _{str(i).rjust(5, '0')}.pdb")
        print(f"_{str(i).rjust(5, '0')}.pdb")
        os.system(f'rm {c_pdb}')

os.system('/Users/rodrigo/software/anaconda3/bin/pdb_mkensemble _*pdb > movie.pdb')


