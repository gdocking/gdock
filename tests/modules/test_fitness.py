import unittest
from modules.fitness import run_dcomplex
from utils.files import get_full_path

data_folder = get_full_path('tests', 'test_data')


class TestFitness(unittest.TestCase):
    def setUp(self):
        self.test_complex = f'{data_folder}/tidy_transformed.pdb'

    def test_run_dcomplex(self):
        observed_energy = run_dcomplex(self.test_complex)
        expected_energy = -4.7
        self.assertEqual(observed_energy, expected_energy)
