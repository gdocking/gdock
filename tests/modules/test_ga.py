import os
import random
import unittest
import io
from tempfile import NamedTemporaryFile

from gdock.modules.files import get_full_path
from gdock.modules.ga import GeneticAlgorithm
from gdock.modules.initialize import Setup

DATA_FOLDER = get_full_path("tests", "test_data")
ETC_FOLDER = get_full_path("etc")


class TestGeneticAlgorithm(unittest.TestCase):
    def setUp(self):
        pioneer_pdb = f"{DATA_FOLDER}/complex.pdb"
        with open(pioneer_pdb) as fh:
            test_pioneer = "".join(fh.readlines())
        fh.close()

        # This could be in a file, but keep it here for "practical" reasons
        toml_string = "[main]" + os.linesep
        toml_string += "identifier = 'setup'" + os.linesep
        toml_string += "number_of_processors = 1" + os.linesep
        toml_string += "random_seed = 42" + os.linesep
        toml_string += "[ga]" + os.linesep
        toml_string += "population_size = 1" + os.linesep
        toml_string += "max_number_of_generations = 1" + os.linesep
        toml_string += "[restraints]" + os.linesep
        toml_string += "A = [39,40,41]" + os.linesep
        toml_string += "B = [4,5,6]" + os.linesep
        toml_string += "[molecules]" + os.linesep
        toml_string += f"A = '{DATA_FOLDER}/molA.pdb'" + os.linesep
        toml_string += f"B = '{DATA_FOLDER}/molB.pdb'" + os.linesep

        test_param_f = NamedTemporaryFile(delete=False, suffix=".toml")
        test_param_f.write(str.encode(toml_string))
        test_param_f.close()
        test_run_params = Setup(test_param_f.name).initialize()

        self.GeneticAlgorithm = GeneticAlgorithm(test_pioneer, test_run_params)
        self.GeneticAlgorithm.setup()

    def test_setup(self):
        observed_toolbox = self.GeneticAlgorithm.toolbox
        self.assertTrue(observed_toolbox)

    def test_run(self):
        observed_ga_dic = self.GeneticAlgorithm.run()

        self.assertTrue("individual" in observed_ga_dic[1][0])
        self.assertTrue("fitness" in observed_ga_dic[1][0])
        self.assertTrue("clone" in observed_ga_dic[1][0])

        observed_individual = observed_ga_dic[1][0]["individual"]
        observed_fitness = observed_ga_dic[1][0]["fitness"]
        observed_structure = observed_ga_dic[1][0]["structure"]
        observed_clone = observed_ga_dic[1][0]["clone"]

        self.assertEqual(observed_individual, [327, 57, 12, 0, -1, -1])
        self.assertEqual(round(observed_fitness[0], 2), -1.88)
        self.assertTrue("0001_0000.pdb" in observed_structure)
        self.assertIsNone(observed_clone)

    def test_fitness_function(self):
        test_pdb_dic = self.GeneticAlgorithm.pioneer_dic
        test_individual = [0, 0, 0, 0, 0, 0]
        # test dcomplex
        observed_fitness_l = self.GeneticAlgorithm.fitness_function(
            test_pdb_dic, "dcomplex", test_individual
        )
        self.assertEqual(round(observed_fitness_l[0], 2), 1.27)
        # test haddock-score
        observed_fitness_l = self.GeneticAlgorithm.fitness_function(
            test_pdb_dic, "haddock", test_individual
        )
        self.assertEqual(round(observed_fitness_l[0], 2), -8.08)

    def test_generate_individual(self):
        random.seed(42)
        observed_individual = self.GeneticAlgorithm.generate_individual()
        expected_individual = [327, 57, 12, 0, -1, -1]
        self.assertEqual(observed_individual, expected_individual)

    # def test__generate_epoch(self):
    #     pass

    def test__recreate(self):
        input_strct_dic = self.GeneticAlgorithm.pioneer_dic
        individual = [327, 57, 12, -1, -2, -2]
        name = f"{DATA_FOLDER}/temp.pdb"
        exected = f"{DATA_FOLDER}/recreated.pdb"

        self.GeneticAlgorithm._recreate(input_strct_dic, individual, name)

        # check if the structure was generated
        self.assertTrue(os.path.isfile(name))

        self.assertListEqual(list(io.open(name)), list(io.open(exected)))

        os.unlink(name)
