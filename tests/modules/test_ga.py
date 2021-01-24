import unittest
import toml
import random
from tempfile import NamedTemporaryFile
from modules.setup import Setup
from modules.ga import GeneticAlgorithm
from utils.files import get_full_path

data_folder = get_full_path('tests', 'test_data')


class TestGeneticAlgorithm(unittest.TestCase):
    def setUp(self):
        test_pioneer = ''.join(open(f'{data_folder}/tidy_transformed.pdb').readlines())
        test_ga_params = toml.load(f'{data_folder}/ga_test_params.toml')
        random.seed(test_ga_params['parameters']['random_seed'])
        toml_string = f"""
[main]
identifier = 'setup'
number_of_processors = 1

[restraints]
A = [39,40,41]
B = [4,5,6]

[molecules]
A = '{data_folder}/molA.pdb'
B = '{data_folder}/molB.pdb'
"""
        test_param_f = NamedTemporaryFile(delete=False, suffix='.toml')
        test_param_f.write(str.encode(toml_string))
        test_param_f.close()
        test_run_params, _ = Setup(test_param_f.name).initialize()

        self.GeneticAlgorithm = GeneticAlgorithm(test_pioneer, test_run_params, test_ga_params)
        self.GeneticAlgorithm.setup()

    def test_setup(self):
        observed_toolbox = self.GeneticAlgorithm.toolbox
        self.assertTrue(observed_toolbox)

    def test_run(self):
        observed_ga_dic = self.GeneticAlgorithm.run()
        expected_ga_dic = {1: {0: ([327, 57, 12, -1, -2, -2], [-4.7])}}
        self.assertDictEqual(observed_ga_dic, expected_ga_dic)

    def test_fitness_function(self):
        test_pdb_dic = self.GeneticAlgorithm.pioneer_dic
        test_individual = [327, 57, 12, -1, -2, -2]
        observed_fitness_l = self.GeneticAlgorithm.fitness_function(test_pdb_dic, test_individual)
        expected_fitness_l = [-4.7]
        self.assertEqual(observed_fitness_l, expected_fitness_l)

    def test_generate_individual(self):
        observed_individual = self.GeneticAlgorithm.generate_individual()
        expected_individual = [327, 57, 12, -1, -2, -2]
        self.assertEqual(observed_individual, expected_individual)
