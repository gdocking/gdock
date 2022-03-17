import os
import pathlib
import shutil
import sys
import unittest
from pathlib import Path

import toml

sys.path.append(str(Path(Path(__file__).parent.parent.parent, "src")))

from gdock.modules.analysis import Analysis
from gdock.modules.files import get_full_path

DATA_FOLDER = get_full_path("tests", "test_data")
ETC_FOLDER = get_full_path("etc")


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
                    "structure": f"{DATA_FOLDER}/0001_0000.pdb",
                    "clone": None,
                    "ranking": 1,
                    "score": 25.06220902040817,
                    "energy": 10.772353,
                }
            }
        }

        structures_path = pathlib.Path(DATA_FOLDER) / "clustering"
        structure_list = structures_path.glob("*pdb")
        self.structure_list = [str(structure) for structure in structure_list]
        self.structure_list.sort()
        self.native = f"{DATA_FOLDER}/1a2k_complex_bound.pdb"
        self.nproc = 1

        params = toml.load(Path(ETC_FOLDER, "params.toml"))
        params["folder"] = str(self.root_folder)
        params["native"] = self.native
        self.Analysis = Analysis(self.result_dic, params)

    def test_get_structures(self):
        observed_l = self.Analysis.get_structures(self.result_dic)

        first_pdb = pathlib.Path(observed_l[0])
        self.assertEqual(first_pdb.stem, "0001_0000")

    def test_cluster(self):

        self.Analysis.structure_list = self.structure_list
        self.Analysis.clust_cutoff = 0.4
        self.Analysis.clust_min_size = 4

        self.Analysis.cluster()

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
        self.Analysis.evaluate_irmsd()

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

        expected_header = "gen\tind\tranking\tfitness\tirmsd\tcluster_id\tinternal_cluster_ranking\tstructure_path\n"

        self.assertEqual(header, expected_header)

    def test_generate_plots(self):
        test_dat = Path(DATA_FOLDER, "test.dat")
        self.Analysis.generate_plots(test_dat)

        kde_plot = Path(DATA_FOLDER, "kde.png")
        ridge_plot = Path(DATA_FOLDER, "ridge.png")

        self.assertTrue(kde_plot.exists())
        self.assertTrue(ridge_plot.exists())

        os.unlink(kde_plot)
        os.unlink(ridge_plot)

        kde_plot = Path(DATA_FOLDER, "kde.svg")
        ridge_plot = Path(DATA_FOLDER, "ridge.svg")

        self.assertTrue(kde_plot.exists())
        self.assertTrue(ridge_plot.exists())

        os.unlink(kde_plot)
        os.unlink(ridge_plot)

    def tearDown(self):
        shutil.rmtree(self.root_folder)
