import numpy as np
from scipy.spatial.transform import Rotation
from utils.functions import tidy  # , write_coords, add_dummy
import logging
import warnings
warnings.filterwarnings('ignore', '.*Optimal rotation is not unique.*')

ga_log = logging.getLogger('ga_log')


class Geometry:

    def __init__(self, input_data, restraint):
        """Initialize the geometry class."""
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

        mat, _ = Rotation.align_vectors(a, b)

        rot_l_c = mat.apply(l_c)
        rot_l_rest_c = mat.apply(l_rest_c)

        rot_l_center = rot_l_c.mean(axis=0)
        rot_l_rest_center = rot_l_rest_c.mean(axis=0)

        # Move them apart
        a = np.array([r_center, r_rest_center_c])
        b = np.array([rot_l_center, rot_l_rest_center])

        _, b_i = np.add(a, b + 5)  # +5 to keep them further apart

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
        rotation = Rotation.from_rotvec(rotation_vector)

        l_c = rotation.apply(rot_l_c)
        l_rest_c = rotation.apply(rot_l_rest_c)

        l_c += c
        l_rest_c += c

        self.ligand_coord = l_c
        self.receptor_coord = r_c

    def apply_transformation(self):
        """Apply transformations and place partners facing each other."""
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
    def translate(coords, point):
        """Translation function."""
        trans_coords = coords + point
        return trans_coords

    @staticmethod
    def rotate(coords, rotation):
        """Rotation function."""
        center = coords.mean(axis=0)
        coords -= center
        rot = Rotation.from_euler('zyx', rotation)
        r = np.array([rot.apply(e) for e in coords])
        r += center
        rotated_coords = r
        return rotated_coords
