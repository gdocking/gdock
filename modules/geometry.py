import numpy as np
from scipy.spatial.transform import Rotation as R
# from utils.functions import add_dummy, write_coords
import logging
ga_log = logging.getLogger('ga_log')

class Geometry:

    def __init__(self, input, restraint):
        self.receptor_coord = input.coords['A']
        self.ligand_coord = input.coords['B']
        self.receptor_rest_coord = restraint.coords['A']
        self.ligand_rest_coord = restraint.coords['B']

    def initial_position(self):
        # calculate the geometric center of the molecule and of the restraints
        ga_log.info('Positioning molecules in starting conformation')
        r_c = np.array(self.receptor_coord)
        l_c = np.array(self.ligand_coord)

        r_rest_c = np.array(self.receptor_rest_coord)
        l_rest_c = np.array(self.ligand_rest_coord)

        r_center = r_c.mean(axis=0)
        l_center = l_c.mean(axis=0)

        # Move both to origin
        r_c -= r_center
        l_c -= l_center

        r_rest_c -= r_center
        l_rest_c -= l_center

        # get center of interface and molecule
        r_center = r_c.mean(axis=0)
        l_center = l_c.mean(axis=0)

        r_rest_center_c = r_rest_c.mean(axis=0)
        l_rest_center_c = l_rest_c.mean(axis=0)

        # write_coords('/Users/rodrigo/repos/gadock/dev/input/target-unbound_A.pdb',
        #              '/Users/rodrigo/repos/gadock/dev/input/A-ori.pdb',
        #              r_c)
        #
        # write_coords('/Users/rodrigo/repos/gadock/dev/input/target-unbound_B.pdb',
        #              '/Users/rodrigo/repos/gadock/dev/input/B-ori.pdb',
        #              l_c)
        #
        # add_dummy('/Users/rodrigo/repos/gadock/dev/input/A-ori.pdb',
        #           '/Users/rodrigo/repos/gadock/dev/input/A-ori-dummy.pdb',
        #           (r_center,r_rest_center_c))
        #
        # add_dummy('/Users/rodrigo/repos/gadock/dev/input/B-ori.pdb',
        #           '/Users/rodrigo/repos/gadock/dev/input/B-ori-dummy.pdb',
        #           (l_center,l_rest_center_c))

        #====
        # Align the vectors between molecule center and interface

        a = np.array([r_center, r_rest_center_c])
        b = np.array([l_center, l_rest_center_c])

        mat, rmsd = R.align_vectors(a, b)

        rot_l_c = mat.apply(l_c)
        rot_l_rest_c = mat.apply(l_rest_c)

        rot_l_center = rot_l_c.mean(axis=0)
        rot_l_rest_center = rot_l_rest_c.mean(axis=0)

        # write_coords('/Users/rodrigo/repos/gadock/dev/input/B-ori.pdb',
        #              '/Users/rodrigo/repos/gadock/dev/input/B-rot.pdb',
        #              rot_l_c)
        #
        # add_dummy('/Users/rodrigo/repos/gadock/dev/input/B-rot.pdb',
        #           '/Users/rodrigo/repos/gadock/dev/input/B-rot-dummy.pdb',
        #           (rot_l_center, rot_l_rest_center))

        # Move them apart
        a = np.array([r_center, r_rest_center_c])
        b = np.array([rot_l_center, rot_l_rest_center])

        a_i, b_i = np.add(a, b)

        rot_l_c += b_i
        rot_l_rest_c += b_i

        # write_coords('/Users/rodrigo/repos/gadock/dev/input/B-rot.pdb',
        #              '/Users/rodrigo/repos/gadock/dev/input/B-trans.pdb',
        #              rot_l_c)

        # Rotate so that they face each other

        # Move to origin and rotate again
        c = rot_l_c.mean(axis=0)
        rot_l_c -= c
        rot_l_rest_c -= c

        rotation_radians = np.radians(180)
        rotation_axis = np.array([0, 0, 1])
        rotation_vector = rotation_radians * rotation_axis
        rotation = R.from_rotvec(rotation_vector)

        l_c = rotation.apply(rot_l_c)
        l_rest_c = rotation.apply(rot_l_rest_c)

        l_c += c
        l_rest_c += c

        l_center = l_c.mean(axis=0)
        l_rest_center = l_rest_c.mean(axis=0)

        # write_coords('/Users/rodrigo/repos/gadock/dev/input/B-rot.pdb',
        #              '/Users/rodrigo/repos/gadock/dev/input/B-maybe.pdb',
        #              l_c)
        #
        # add_dummy('/Users/rodrigo/repos/gadock/dev/input/B-maybe.pdb',
        #           '/Users/rodrigo/repos/gadock/dev/input/B-maybe-dummy.pdb',
        #           (l_center, l_rest_center))


    @staticmethod
    def rotate_molecule(mol, rotation_mat):
        pass

    @staticmethod
    def calc_center(coords):
        coord_array = np.array(coords)
        center = coord_array.mean(axis=0)
        return center

    @staticmethod
    def rotation_matrix_from_vectors(vec1, vec2):
        """ Find the rotation matrix that aligns vec1 to vec2
        :param vec1: A 3d "source" vector
        :param vec2: A 3d "destination" vector
        :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
        """
        a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
        v = np.cross(a, b)
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
        return rotation_matrix
