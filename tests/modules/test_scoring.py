import configparser
import os
import unittest
from pathlib import Path

import toml
from gdock.modules.files import get_full_path
from gdock.modules.scoring import Scoring

etc_folder = get_full_path("etc")
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, "gdock.ini"), encoding="utf-8")
dcomplex_exe = ini.get("third_party", "dcomplex_exe")

DATA_FOLDER = get_full_path("tests", "test_data")
ETC_FOLDER = get_full_path("etc")


class TestScoring(unittest.TestCase):
    def setUp(self):
        data_dic = {
            1: {
                0: {
                    "individual": [327, 57, 12, -1, -2, -2],
                    "fitness": (0.33,),
                    "structure": f"{DATA_FOLDER}/0001_0000.pdb",
                    "clone": None,
                }
            }
        }
        params = toml.load(Path(ETC_FOLDER, "params.toml"))
        params["main"]["number_of_processors"] = 1
        self.Scoring = Scoring(data_dic, params)

    def test_score(self):
        observed_scored_dic = self.Scoring.score()

        self.assertTrue("ranking" in observed_scored_dic[1][0])
        self.assertTrue("score" in observed_scored_dic[1][0])
        self.assertTrue("energy" in observed_scored_dic[1][0])

        observed_ranking = observed_scored_dic[1][0]["ranking"]
        observed_score = observed_scored_dic[1][0]["score"]
        observed_energy = observed_scored_dic[1][0]["energy"]

        self.assertEqual(observed_ranking, 1)
        self.assertEqual(round(observed_score, 2), 32.64)
        self.assertEqual(round(observed_energy, 2), 10.77)

    # def test_calculate_energy(self):
    #     # this is tested by test_score
    #     pass

    # def test_run_dcomplex(self):
    #     # this is tested by test_score
    #     pass
