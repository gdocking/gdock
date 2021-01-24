import unittest
import tempfile
import shutil
from utils.functions import tidy, format_coords, du
from utils.files import get_full_path

data_folder = get_full_path('tests', 'test_data')


class TestFunctions(unittest.TestCase):
    def setUp(self):
        with open(f'{data_folder}/molA.pdb', 'r') as test_fh:
            self.pdb_str = ''.join(test_fh.readlines())
        test_fh.close()

        self.coords = [10, -11, 45.2]

        self.temp_dir = tempfile.TemporaryDirectory()
        shutil.copy(f'{data_folder}/molA.pdb', self.temp_dir.name)

    def test_tidy(self):
        observed_tidy_str = tidy(self.pdb_str)
        expected_tidy_str = 'ATOM      1  CA  CYS A   1       2.945 -15.164  18.823  1.00 13.76      A   19  \n' \
                            'ATOM      2  CA  GLY A   2       6.035 -13.796  20.667  1.00 14.70      A   19  \n' \
                            'ATOM      3  CA  VAL A   3       8.279 -16.791  19.929  1.00 15.89      A   19  \n' \
                            'ATOM      4  CA  PRO A   4      10.599 -16.062  16.976  1.00 11.38      A   19  \n' \
                            'ATOM      5  CA  ALA A   5      11.893 -18.739  14.665  1.00 15.49      A   19  \n' \
                            'TER       6      ALA A   5                                                      \n' \
                            'END                                                                             \n'
        self.assertEqual(observed_tidy_str, expected_tidy_str)

    def test_format_coords(self):
        observed_formated_coords = format_coords(self.coords)
        expected_formated_coords = (' 10.000', '-11.000', ' 45.200')
        self.assertEqual(observed_formated_coords, expected_formated_coords)

    def test_random_quote(self):
        # does it make sense to test this one?
        pass

    def test_du(self):
        observed_du = du(self.temp_dir.name)
        expected_du = 567
        self.assertAlmostEqual(observed_du, expected_du)
