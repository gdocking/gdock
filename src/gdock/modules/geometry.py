import logging
import warnings

import numpy as np
from scipy.spatial.transform import Rotation

from gdock.modules.functions import tidy  # , write_coords, add_dummy

warnings.filterwarnings("ignore", ".*Optimal rotation is not unique.*")

ga_log = logging.getLogger("ga_log")


def translate(coords, point):
    """Translation function."""
    trans_coords = coords + point
    return trans_coords


def rotate(coords, rotation):
    """Rotation function."""
    center = coords.mean(axis=0)
    coords -= center
    rot = Rotation.from_euler("zyx", rotation)
    r = np.array([rot.apply(e) for e in coords])
    r += center
    rotated_coords = r
    return rotated_coords


class Complex:
    def __init__(
        self,
        identifier,
        coords_i,
        coords_j,
        restraint_i,
        restraint_j,
        pdb_i,
        pdb_j,
        restraint_i_resnums,
        restraint_j_resnums,
    ):
        self.id = identifier
        self.coords_i = coords_i
        self.coords_j = coords_j
        self.restraint_i = restraint_i
        self.restraint_j = restraint_j
        self.pdb_i = pdb_i
        self.pdb_j = pdb_j
        self.restraint_i_resnums = restraint_i_resnums
        self.restraint_j_resnums = restraint_j_resnums
        # self.receptor = ""
        # self.ligand = ""
        self.structure = None

    def position(self):
        """Position the molecules in the initial position."""
        # ga_log.info("Positioning molecules in starting conformation")

        r_c = np.array(self.coords_i)
        l_c = np.array(self.coords_j)

        r_rest_c = np.array(self.restraint_i)
        l_rest_c = np.array(self.restraint_j)

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

        self.coords_j = l_c
        self.coords_i = r_c

        # ga_log.info("Applying transformations for initial position")

        # ga_log.debug("Applying transformation to the receptor")
        receptor = ""
        for coord, line in zip(self.coords_i, self.pdb_i):
            if line.startswith("ATOM"):
                new_x = f"{coord[0]:.3f}".rjust(7, " ")
                new_y = f"{coord[1]:.3f}".rjust(7, " ")
                new_z = f"{coord[2]:.3f}".rjust(7, " ")
                new_line = f"{line[:30]} {new_x} {new_y} {new_z} {line[55:]}"
                receptor += new_line

        # ga_log.debug("Applying transformation to the ligand")
        ligand = ""
        for coord, line in zip(self.coords_j, self.pdb_j):
            if line.startswith("ATOM"):
                new_x = f"{coord[0]:.3f}".rjust(7, " ")
                new_y = f"{coord[1]:.3f}".rjust(7, " ")
                new_z = f"{coord[2]:.3f}".rjust(7, " ")
                new_line = f"{line[:30]} {new_x} {new_y} {new_z} {line[55:]}"
                ligand += new_line

        tidy_complex = tidy(receptor + ligand)

        self.structure = tidy_complex

    def __repr__(self):
        return f"Complex {self.id}"


class Geometry:
    def __init__(self, input_data, restraints):
        """Initialize the geometry class."""
        self.input = input_data
        self.restraints = restraints
        self.complex_dic = {}

    def prepare_initial_complexes(self):
        """Put the complexes in their initial positions."""
        # input_data.coord might have has multiple models, combine them all
        complex_l = []
        complex_counter = 1
        for model_a in self.input.coords["A"]:
            for model_b in self.input.coords["B"]:
                complex = Complex(
                    complex_counter,
                    self.input.coords["A"][model_a],
                    self.input.coords["B"][model_b],
                    self.restraints.coords["A"][model_a],
                    self.restraints.coords["B"][model_b],
                    self.input.raw_pdb["A"][model_a],
                    self.input.raw_pdb["B"][model_b],
                    self.restraints.resnums["A"],
                    self.restraints.resnums["B"],
                )
                complex.position()

                complex_l.append(complex)

                complex_counter += 1

        return complex_l
