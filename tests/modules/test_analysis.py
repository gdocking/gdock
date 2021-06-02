import unittest
# import configparser
import pathlib
from utils.files import get_full_path
from modules.analysis import Analysis

# etc_folder = get_full_path('etc')
# ini = configparser.ConfigParser(os.environ)
# ini.read(os.path.join(etc_folder, 'gdock.ini'), encoding='utf-8')
# profit_exe = ini.get('third_party', 'profit_exe')

data_folder = get_full_path('tests', 'test_data')


class TestAnalysis(unittest.TestCase):

    def setUp(self):
        folder = pathlib.Path.cwd()
        (folder / 'analysis').mkdir(parents=True, exist_ok=True)
        self.result_dic = {1: {0: {'individual': [327, 57, 12, -1, -2, -2],
                                   'fitness': (0.3333333333333333,),
                                   'structure': f'{data_folder}/0001_0000.pdb',
                                   'clone': None,
                                   'ranking': 1,
                                   'score': 25.06220902040817,
                                   'energy': 10.772353}}}
        run_params = {'np': 1,
                      'native': f'{data_folder}/1a2k_complex_bound.pdb',
                      'folder': str(folder)}

        self.Analysis = Analysis(self.result_dic, run_params)

    def test_get_structures(self):
        observed_l = self.Analysis.get_structures(self.result_dic)

        first_pdb = pathlib.Path(observed_l[0])
        self.assertEqual(first_pdb.stem, '0001_0000')

    # def test_cluster(self):
    #     pass

    # def test_evaluate(self):
    #     pass

    def test_output(self):
        self.Analysis.output()
        output_f = pathlib.Path(self.Analysis.analysis_path) / 'gdock.dat'
        self.assertTrue(output_f.exists())

        with output_f.open() as out_fh:
            header = out_fh.readlines()[0]

        expected_header_elements = ['gen', 'ind', 'ranking', 'score',
                                    'fitness', 'energy', 'irmsd',
                                    'cluster_id', 'internal_cluster_ranking',
                                    'structure_path']

        for expected_element in expected_header_elements:
            self.assertTrue(expected_element in header)
