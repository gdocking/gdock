import numpy as np
from pyquaternion import Quaternion
from utils.functions import timeit, format_coords


class PDB:

    def __init__(self):
        self.dic = {}

    def load(self, pdb_f):
        #
        with open(pdb_f, 'r') as fh:
            for l in fh.readlines():
                if l.startswith('ATOM'):
                    chain = l[21]
                    x = float(l[31:38])
                    y = float(l[39:46])
                    z = float(l[47:54])
                    if chain not in self.dic:
                        self.dic[chain] = {'coord': [], 'raw': []}
                    self.dic[chain]['coord'].append((x, y, z))
                    self.dic[chain]['raw'].append(l)
        fh.close()
        return self.dic

    def rotate(self, target_chain, rotation):
        q = Quaternion(rotation)
        c = np.array(self.dic[target_chain]['coord'])
        center = c.mean(axis=0)
        c -= center
        r = np.array([q.rotate(e) for e in c])
        r += center
        # self.dic[target_chain]['coord'] = list(r)
        return True

    @timeit
    def randomize_rotation(self, target_chain):
        """

        :param target_chain: str
        :return: array of randomly rotated coordinates
        """
        q = Quaternion.random()
        # c = np.array([j[0] for j in self.dic[target_chain]])
        c = np.array(self.dic[target_chain]['coord'])
        center = c.mean(axis=0)
        c -= center

        r = np.array([q.rotate(e) for e in c])
        r += center

        self.dic[target_chain]['coord'] = list(r)

        return list(r)

    def output(self, output_fname):
        with open(output_fname, 'w') as out_fh:
            for chain in self.dic:
                for coord, line in zip(self.dic[chain]['coord'], self.dic[chain]['raw']):
                    new_x, new_y, new_z = format_coords(coord)
                    new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}'
                    out_fh.write(new_line)
        out_fh.close()
        return True
