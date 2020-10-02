# import numpy as np
# from dual_quaternions import DualQuaternion
# from pyquaternion import Quaternion
# from utils.functions import timeit, format_coords
import logging

ga_log = logging.getLogger('ga_log')


class PDB:

    def __init__(self):
        """

        """
        self.coords = {}
        self.raw_pdb = {}

    def load(self, pdb_f):
        ga_log.debug(f'Loading {pdb_f}')
        with open(pdb_f, 'r') as fh:
            for line in fh.readlines():
                if line.startswith('ATOM'):
                    chain = line[21]
                    if chain not in self.coords:
                        self.coords[chain] = []
                    if chain not in self.raw_pdb:
                        self.raw_pdb[chain] = []
                    x = float(line[31:38])
                    y = float(line[39:46])
                    z = float(line[47:54])
                    self.coords[chain].append((x, y, z))
                    self.raw_pdb[chain].append(line)
        fh.close()

    # def rotate(self, target_chain, rotation):
    #     q = Quaternion(rotation)
    #     c = np.array(self.dic[target_chain]['coord'])
    #     center = c.mean(axis=0)
    #     c -= center
    #     r = np.array([q.rotate(e) for e in c])
    #     r += center
    #     # self.dic[target_chain]['coord'] = list(r)
    #     return True
    #
    # def randomize(self, target_chain):
    #     dq = DualQuaternion(Quaternion.random(), Quaternion.random())
    #     c = np.array(self.dic[target_chain]['coord'])
    #     center = c.mean(axis=0)
    #     c -= center
    #
    #     r = np.array([dq.transform_point(e) for e in c])
    #     r += center
    #
    #     self.dic[target_chain]['coord'] = list(r)
    #
    #     return list(r)
    # def centralize(self, target_chain):
    #     c_a = np.array(self.dic['A']['coord'])
    #     c_b = np.array(self.dic[target_chain]['coord'])
    #
    #     center_a = c_a.mean(axis=0)
    #     center_b = c_b.mean(axis=0)
    #     c_b -= [7.657,  -0.102,  17.979]
    #
    #     self.dic[target_chain]['coord'] = list(c_b)
    #
    #     return list(c_b)
    # @timeit
    # def randomize_rotation(self, target_chain):
    #     """
    #
    #     :param target_chain: str
    #     :return: array of randomly rotated coordinates
    #     """
    #     q = Quaternion.random()
    #     # c = np.array([j[0] for j in self.dic[target_chain]])
    #     c = np.array(self.dic[target_chain]['coord'])
    #     center = c.mean(axis=0)
    #     c -= center
    #
    #     r = np.array([q.rotate(e) for e in c])
    #     r += center
    #
    #     self.dic[target_chain]['coord'] = list(r)
    #
    #     return list(r)
    # def output(self, output_fname):
    #     with open(output_fname, 'w') as out_fh:
    #         for chain in self.dic:
    #             for coord, line in zip(self.dic[chain]['coord'], self.dic[chain]['raw']):
    #                 new_x, new_y, new_z = format_coords(coord)
    #                 new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}'
    #                 out_fh.write(new_line)
    #     out_fh.close()
    #     return True


class Restraint:
    def __init__(self, raw_pdb):
        """

        :param raw_pdb:
        """
        self.raw_pdb = raw_pdb
        self.coords = {}

    def load(self, restraint, identifier):
        ga_log.debug(f'Loading restraints {identifier}, {restraint}')
        if identifier not in self.coords:
            self.coords[identifier] = []

        for line in self.raw_pdb[identifier]:
            resnum = int(line[22:26])
            if resnum in restraint:
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                self.coords[identifier].append((x, y, z))
