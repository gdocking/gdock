import unittest
from modules.structure import PDB, Restraint
from utils.files import get_full_path

data_folder = get_full_path('test', 'test_data')


class TestPDB(unittest.TestCase):
    def setUp(self):
        self.PDB = PDB()

    def test_load(self):
        pass

class TestRestraint(unittest.TestCase):
    def setUp(self):
        raw_pdb = {'A': ['ATOM      1  CA  CYS A   1       2.945 -15.164  18.823  1.00 13.76      A   19  \n',
                         'ATOM      2  CA  GLY A   2       6.035 -13.796  20.667  1.00 14.70      A   19  \n',
                         'ATOM      3  CA  VAL A   3       8.279 -16.791  19.929  1.00 15.89      A   19  \n',
                         'ATOM      4  CA  PRO A   4      10.599 -16.062  16.976  1.00 11.38      A   19  \n',
                         'ATOM      5  CA  ALA A   5      11.893 -18.739  14.665  1.00 15.49      A   19  \n'],
                   'B': ['ATOM      1  CA  CYS B   1       2.945 -15.164  18.823  1.00 13.76      A   19  \n',
                         'ATOM      2  CA  GLY B   2       6.035 -13.796  20.667  1.00 14.70      A   19  \n',
                         'ATOM      3  CA  VAL B   3       8.279 -16.791  19.929  1.00 15.89      A   19  \n',
                         'ATOM      4  CA  PRO B   4      10.599 -16.062  16.976  1.00 11.38      A   19  \n',
                         'ATOM      5  CA  ALA B   5      11.893 -18.739  14.665  1.00 15.49      A   19  \n']}
        self.Restraint = Restraint(raw_pdb=raw_pdb)


    def test_load(self):
        expected = {'A': [(6.035, -13.796, 20.667),
                          (8.279, -16.791, 19.929),
                          (10.599, -16.062, 16.976)]}
        self.Restraint.load(restraint=[2,3,4], identifier='A')
        self.assertEqual(self.Restraint.coords, expected)
