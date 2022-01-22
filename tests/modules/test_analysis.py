import pathlib
import shutil
import unittest
import sys
from pathlib import Path

sys.path.append(str(Path(Path(__file__).parent.parent.parent, "src")))

from gdock.modules.analysis import Analysis
from gdock.modules.files import get_full_path

data_folder = get_full_path("tests", "test_data")


class TestAnalysis(unittest.TestCase):
    def setUp(self):
        self.root_folder = pathlib.Path.cwd() / "tests" / "dummy_f"
        self.analysis_folder = self.root_folder / "analysis"
        self.analysis_folder.mkdir(parents=True, exist_ok=True)

        self.result_dic = {
            1: {
                0: {
                    "individual": [327, 57, 12, -1, -2, -2],
                    "fitness": (0.3333333333333333,),
                    "structure": f"{data_folder}/0001_0000.pdb",
                    "clone": None,
                    "ranking": 1,
                    "score": 25.06220902040817,
                    "energy": 10.772353,
                }
            }
        }

        structures_path = pathlib.Path(data_folder) / "clustering"
        structure_list = structures_path.glob("*pdb")
        self.structure_list = [str(structure) for structure in structure_list]
        self.structure_list.sort()
        self.native = f"{data_folder}/1a2k_complex_bound.pdb"
        self.nproc = 1

        run_params = {
            "np": self.nproc,
            "native": self.native,
            "folder": str(self.root_folder),
        }

        self.Analysis = Analysis(self.result_dic, run_params)

    def test_get_structures(self):
        observed_l = self.Analysis.get_structures(self.result_dic)

        first_pdb = pathlib.Path(observed_l[0])
        self.assertEqual(first_pdb.stem, "0001_0000")

    def test_cluster(self):

        self.Analysis.structure_list = self.structure_list

        self.Analysis.cluster(cutoff=0.4, min_size=4)

        for structure in self.structure_list:
            contact_f = pathlib.Path(structure.replace(".pdb", ".contacts"))
            contact_f.unlink(missing_ok=True)

        cluster_out = self.analysis_folder / "cluster.out"
        fcc_matrix = self.analysis_folder / "fcc.matrix"

        self.assertTrue(cluster_out.exists())
        self.assertTrue(fcc_matrix.exists())

        with cluster_out.open() as fh:
            observed_cluster_l = fh.readlines()
        fh.close()

        observed_cluster_line_1 = " ".join(observed_cluster_l[0].split())
        observed_cluster_line_2 = " ".join(observed_cluster_l[1].split())

        self.assertEqual(observed_cluster_line_1, "Cluster 1 -> 9 27 28 31 32")
        self.assertEqual(observed_cluster_line_2, "Cluster 2 -> 8 14 25 26 29")

    def test_evaluate(self):
        self.Analysis.structure_list = self.structure_list[:3]
        self.Analysis.evaluate()

        observed_irmsd_dic = self.Analysis.irmsd_dic
        expected_irmsd_dic = {
            "0001_0000": 14.196,
            "0001_0001": 17.489,
            "0001_0002": 15.42,
        }
        self.assertDictEqual(observed_irmsd_dic, expected_irmsd_dic)

    def test_output(self):
        self.Analysis.output()
        output_f = self.analysis_folder / "gdock.dat"
        self.assertTrue(output_f.exists())

        with output_f.open() as out_fh:
            header = out_fh.readlines()[0]

        expected_header_elements = [
            "gen",
            "ind",
            "ranking",
            "score",
            "fitness",
            "energy",
            "irmsd",
            "cluster_id",
            "internal_cluster_ranking",
            "structure_path",
        ]

        for expected_element in expected_header_elements:
            self.assertTrue(expected_element in header)

    def tearDown(self):
        shutil.rmtree(self.root_folder)
