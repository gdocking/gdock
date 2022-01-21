import os
import unittest

import numpy as np

from modules.geometry import Geometry
from modules.structure import PDB, Restraint
from utils.files import get_full_path

data_folder = get_full_path("tests", "test_data")


class TestGeometry(unittest.TestCase):
    def setUp(self):
        test_input_data = PDB()
        test_input_data.load(f"{data_folder}/molA.pdb")
        test_input_data.load(f"{data_folder}/molB.pdb")

        test_restraints = Restraint(test_input_data.raw_pdb)
        test_restraints.load([1, 2, 3], "A")
        test_restraints.load([2, 3, 4], "B")

        self.Geometry = Geometry(test_input_data, test_restraints)

    def test_calc_initial_position(self):
        self.Geometry.calc_initial_position()
        expected_ligand_coord = np.array(
            [
                [4.41488596, 6.18960229, 2.7898511],
                [5.38935079, 4.86187319, 6.26942201],
                [3.13271966, 7.23145808, 8.22975065],
                [-0.11367046, 5.48686873, 9.2551331],
                [-3.37010188, 7.31580583, 9.73693127],
            ]
        )
        expected_receptor_coord = np.array(
            [
                [-5.0052, 0.9464, 0.611],
                [-1.9152, 2.3144, 2.455],
                [0.3288, -0.6806, 1.717],
                [2.6488, 0.0484, -1.236],
                [3.9428, -2.6286, -3.547],
            ]
        )

        observed_ligand_coord = self.Geometry.ligand_coord
        observed_receptor_coord = self.Geometry.receptor_coord

        self.assertTrue(np.allclose(observed_ligand_coord, expected_ligand_coord))
        self.assertTrue(np.allclose(observed_receptor_coord, expected_receptor_coord))

    def test_apply_transformation(self):
        observed_tidy_complex = self.Geometry.apply_transformation()
        observed_tidy_complex = observed_tidy_complex.split(os.linesep)

        with open(f"{data_folder}/tidy_transformed.pdb", "r") as test_tidy_fh:
            expected_tidy_complex = "".join(test_tidy_fh.readlines())
        test_tidy_fh.close()

        expected_tidy_complex = expected_tidy_complex.split(os.linesep)
        self.assertEqual(len(observed_tidy_complex), len(expected_tidy_complex))
        for observed_line, expected_line in zip(
            observed_tidy_complex, expected_tidy_complex
        ):
            self.assertEqual(observed_line, expected_line)

    def test_translate(self):
        coords = np.array([[34.082, -15.413, 35.895], [34.912, -15.074, 34.747]])
        point = [-1, +2, -3]
        observed_trans_coords = self.Geometry.translate(coords, point)
        expected_trans_coords = np.array(
            [[33.082, -13.413, 32.895], [33.912, -13.074, 31.747]]
        )
        self.assertTrue(np.allclose(observed_trans_coords, expected_trans_coords))

    def test_rotate(self):
        coords = np.array([[34.082, -15.413, 35.895], [34.912, -15.074, 34.747]])
        rotation = [327, 57, 12]
        observed_rot_coords = self.Geometry.rotate(coords, rotation)
        expected_rot_coords = np.array(
            [
                [34.42920652, -15.11616826, 36.03487809],
                [34.56479348, -15.37083174, 34.60712191],
            ]
        )

        self.assertTrue(np.allclose(observed_rot_coords, expected_rot_coords))
