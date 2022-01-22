import unittest

# from modules.fitness import run_dcomplex
from gdock.modules.files import get_full_path

data_folder = get_full_path("tests", "test_data")


class TestFitness(unittest.TestCase):
    def setUp(self):
        self.test_complex = f"{data_folder}/tidy_transformed.pdb"

    def test_calc_satisfaction(self):
        pass
