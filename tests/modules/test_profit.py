import configparser
import math
import os
import tempfile
import unittest

from modules.profit import Profit
from utils.files import get_full_path

etc_folder = get_full_path("etc")
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, "gdock.ini"), encoding="utf-8")
profit_exe = ini.get("third_party", "profit_exe")

data_folder = get_full_path("tests", "test_data")


class TestProfit(unittest.TestCase):
    def setUp(self):
        mobi = [f"{data_folder}/0001_0000.pdb"]
        ref = f"{data_folder}/1a2k_complex_bound.pdb"
        nproc = 1

        self.Profit = Profit(ref, mobi, nproc)

    def test_calc_irmsd(self):
        observed_irmsd_dic = self.Profit.calc_irmsd()

        self.assertTrue("0001_0000" in observed_irmsd_dic)
        self.assertEqual(round(observed_irmsd_dic["0001_0000"], 3), 6.138)

    def test__write_script(self):
        mobi = "some_pdb.pdb"
        observed_script_f = self.Profit._write_script(mobi)
        observed_script_str = ""
        with open(observed_script_f.name, "r") as fh:
            for line in fh.readlines():
                observed_script_str += line

        expected_script_str = f"ref {self.Profit.reference}" + os.linesep
        expected_script_str += f"mobi {mobi}" + os.linesep
        expected_script_str += "atoms C,CA,N,O" + os.linesep
        expected_script_str += self.Profit.izone
        expected_script_str += "fit" + os.linesep
        expected_script_str += "quit"

        self.assertEqual(observed_script_str, expected_script_str)

    def test_execute(self):
        mobi_pdb = f"{data_folder}/0001_0000.pdb"
        ref_pdb = f"{data_folder}/1a2k_complex_bound.pdb"

        script_f = f"ref {ref_pdb}" + os.linesep
        script_f += f"mobi {mobi_pdb}" + os.linesep
        script_f += "atoms C,CA,N,O" + os.linesep
        script_f += self.Profit.izone
        script_f += "fit" + os.linesep
        script_f += "quit"

        temp_script = tempfile.NamedTemporaryFile(delete=False, suffix=".txt")
        temp_script.write(str.encode(script_f))
        temp_script.close()

        value = self.Profit.execute(profit_exe, temp_script.name)
        self.assertEqual(value, 6.138)

        os.unlink(temp_script.name)

        value = self.Profit.execute(profit_exe, "")
        self.assertTrue(math.isnan(value))
