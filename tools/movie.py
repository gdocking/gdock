# read the output file and create a trajectory based on the fitness
import os
from operator import itemgetter

movie_l = []
with open('gadock.dat','r') as fh:
    for l in fh.readlines():
        data = l.split(',')
        gen = data[0]
        ind = data[1]
        fitness = float(data[2])
        pdb = data[3].split('\n')[0]
        movie_l.append((pdb, fitness))


sorted_movie_l = sorted(movie_l, key=itemgetter(1))
sorted_movie_l.reverse()

for i, e in enumerate(sorted_movie_l):
    pdb, fitness = e
    os.system(f"/Users/rodrigo/software/anaconda3/bin/pdb_tidy {pdb} > _{str(i).rjust(5, '0')}.pdb")

os.system('/Users/rodrigo/software/anaconda3/bin/pdb_mkensemble _*pdb > movie.pdb')


