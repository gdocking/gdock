import numpy as np
from scipy.spatial.transform import Rotation as R
from utils.functions import tidy  # , write_coords, add_dummy
import logging
ga_log = logging.getLogger('ga_log')


class Geometry:

    def __init__(self, input_data, restraint):
        """

        :param input_data:
        :param restraint:
        """
        self.receptor_coord = input_data.coords['A']
        self.receptor_pdb = input_data.raw_pdb['A']
        self.receptor_rest_coord = restraint.coords['A']
        self.begin_receptor = ''
        self.ligand_coord = input_data.coords['B']
        self.ligand_pdb = input_data.raw_pdb['B']
        self.ligand_rest_coord = restraint.coords['B']
        self.begin_ligand = ''

    def calc_initial_position(self):
        """Position the molecules in the initial position."""

        # calculate the geometric center of the molecule and of the restraints
        # Q: Maybe use quaternions here as well!

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

        # ====
        # Align the vectors between molecule center and interface

        a = np.array([r_center, r_rest_center_c])
        b = np.array([l_center, l_rest_center_c])

        mat, _ = R.align_vectors(a, b)

        rot_l_c = mat.apply(l_c)
        rot_l_rest_c = mat.apply(l_rest_c)

        rot_l_center = rot_l_c.mean(axis=0)
        rot_l_rest_center = rot_l_rest_c.mean(axis=0)

        # Move them apart
        a = np.array([r_center, r_rest_center_c])
        b = np.array([rot_l_center, rot_l_rest_center])

        _, b_i = np.add(a, b)

        rot_l_c += b_i
        rot_l_rest_c += b_i

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

        self.ligand_coord = l_c
        self.receptor_coord = r_c

    def apply_transformation(self):
        """Apply transformations to put the binding partners in the appropriate places."""

        ga_log.info('Applying transformations for initial position')

        ga_log.debug('Applying transformation to the receptor')
        for coord, line in zip(self.receptor_coord, self.receptor_pdb):
            if line.startswith('ATOM'):
                new_x = f'{coord[0]:.3f}'.rjust(7, ' ')
                new_y = f'{coord[1]:.3f}'.rjust(7, ' ')
                new_z = f'{coord[2]:.3f}'.rjust(7, ' ')
                new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}'
                self.begin_receptor += new_line

        ga_log.debug('Applying transformation to the ligand')
        for coord, line in zip(self.ligand_coord, self.ligand_pdb):
            if line.startswith('ATOM'):
                new_x = f'{coord[0]:.3f}'.rjust(7, ' ')
                new_y = f'{coord[1]:.3f}'.rjust(7, ' ')
                new_z = f'{coord[2]:.3f}'.rjust(7, ' ')
                new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}'
                self.begin_ligand += new_line

        tidy_complex = tidy(self.begin_receptor + self.begin_ligand)

        return tidy_complex

    @staticmethod
    def calc_center(coords):
        """Calculate the geometric center."""
        coord_array = np.array(coords)
        center = coord_array.mean(axis=0)
        return center

    @staticmethod
    def rotation_matrix_from_vectors(vec1, vec2):
        """
        Find the rotation matrix that aligns vec1 to vec2.

        :param vec1: A 3d "source" vector
        :param vec2: A 3d "destination" vector
        :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
        """
        # Copied from https://stackoverflow.com/a/59204638
        a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
        v = np.cross(a, b)
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
        return rotation_matrix
