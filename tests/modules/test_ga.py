import unittest
import toml
import random
import os
from tempfile import NamedTemporaryFile
from modules.setup import Setup
from modules.ga import GeneticAlgorithm
from utils.files import get_full_path

data_folder = get_full_path('tests', 'test_data')


class TestGeneticAlgorithm(unittest.TestCase):
    def setUp(self):
        pioneer_pdb = f'{data_folder}/tidy_transformed.pdb'
        test_pioneer = ''.join(open(pioneer_pdb).readlines())
        test_ga_params = toml.load(f'{data_folder}/ga_test_params.toml')
        random.seed(test_ga_params['parameters']['random_seed'])

        toml_string = "[main]" + os.linesep
        toml_string += "identifier = 'setup'" + os.linesep
        toml_string += "number_of_processors = 1" + os.linesep
        toml_string += "[restraints]" + os.linesep
        toml_string += "A = [39,40,41]" + os.linesep
        toml_string += "B = [4,5,6]" + os.linesep
        toml_string += "[molecules]" + os.linesep
        toml_string += f"A = '{data_folder}/molA.pdb'" + os.linesep
        toml_string += f"B = '{data_folder}/molB.pdb'" + os.linesep

        test_param_f = NamedTemporaryFile(delete=False, suffix='.toml')
        test_param_f.write(str.encode(toml_string))
        test_param_f.close()
        test_run_params, _ = Setup(test_param_f.name).initialize()

        self.GeneticAlgorithm = GeneticAlgorithm(test_pioneer, test_run_params,
                                                 test_ga_params)
        self.GeneticAlgorithm.setup()

    def test_setup(self):
        observed_toolbox = self.GeneticAlgorithm.toolbox
        self.assertTrue(observed_toolbox)

    def test_run(self):
        observed_ga_dic = self.GeneticAlgorithm.run()

        self.assertTrue('individual' in observed_ga_dic[1][0])
        self.assertTrue('fitness' in observed_ga_dic[1][0])
        self.assertTrue('clone' in observed_ga_dic[1][0])

        observed_individual = observed_ga_dic[1][0]['individual']
        observed_fitness = observed_ga_dic[1][0]['fitness']
        observed_structure = observed_ga_dic[1][0]['structure']
        observed_clone = observed_ga_dic[1][0]['clone']

        self.assertEqual(observed_individual, [327, 57, 12, -1, -2, -2])
        self.assertEqual(observed_fitness, (0.0,))
        self.assertTrue('0001_0000.pdb' in observed_structure)
        self.assertIsNone(observed_clone)

    def test_fitness_function(self):
        test_pdb_dic = self.GeneticAlgorithm.pioneer_dic
        test_individual = [327, 57, 12, -1, -2, -2]
        observed_fitness_l = self.GeneticAlgorithm.fitness_function(
            test_pdb_dic, test_individual)
        expected_fitness_l = [0.0]
        self.assertEqual(observed_fitness_l, expected_fitness_l)

    def test_generate_individual(self):
        observed_individual = self.GeneticAlgorithm.generate_individual()
        expected_individual = [327, 57, 12, -1, -2, -2]
        self.assertEqual(observed_individual, expected_individual)

    # def test__generate_epoch(self):
    #     pass

    def test__recreate(self):
        input_strct_dic = self.GeneticAlgorithm.pioneer_dic
        individual = [327, 57, 12, -1, -2, -2]
        name = f'{data_folder}/temp.pdb'

        self.GeneticAlgorithm._recreate(input_strct_dic, individual, name)

        self.assertTrue(os.path.isfile(name))

        expected_file_str = ('ATOM      1  CA  CYS A   1       2.945 -15.164'
                             '  18.823  1.00 13.76      A   19  ' + os.linesep)
        expected_file_str += ('ATOM      2  CA  GLY A   2       6.035 -13.796'
                              '  20.667  1.00 14.70      A   '
                              '19  ' + os.linesep)
        expected_file_str += ('ATOM      3  CA  VAL A   3       8.279 -16.791'
                              '  19.929  1.00 15.89      A   '
                              '19  ' + os.linesep)
        expected_file_str += ('ATOM      4  CA  PRO A   4      10.599 -16.062'
                              '  16.976  1.00 11.38      A   '
                              '19  ' + os.linesep)
        expected_file_str += ('ATOM      5  CA  ALA A   5      11.893 -18.739'
                              '  14.665  1.00 15.49      A   '
                              '19  ' + os.linesep)
        expected_file_str += ('ATOM      2  CA  CYS B   1       2.650 -17.003'
                              '  18.782  1.00 13.76      A   '
                              '19  ' + os.linesep)
        expected_file_str += ('ATOM      3  CA  GLY B   2       5.798 -14.905'
                              '  18.069  1.00 14.70      A   '
                              '19  ' + os.linesep)
        expected_file_str += ('ATOM      4  CA  VAL B   3       8.150 -17.876'
                              '  17.635  1.00 15.89      A   '
                              '19  ' + os.linesep)
        expected_file_str += ('ATOM      5  CA  PRO B   4       8.694 -18.656'
                              '  13.930  1.00 11.38      A   '
                              '19  ' + os.linesep)
        expected_file_str += ('ATOM      6  CA  ALA B   5       9.459 -22.112'
                              '  12.644  1.00 15.49      A   '
                              '19  ' + os.linesep)

        observed_file_str = ''.join(open(name).readlines())
        self.assertEqual(observed_file_str, expected_file_str)
