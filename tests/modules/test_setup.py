import unittest
import pathlib
import shutil
import glob
from modules.setup import Setup
from utils.files import get_full_path

data_folder = get_full_path('tests', 'test_data')


class TestSetup(unittest.TestCase):
    def setUp(self):
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
        with open(f'{data_folder}/setup.toml', 'w') as setup_f:
            setup_f.write(toml_string)
        self.setup_rundir = pathlib.Path(f"{pathlib.Path.cwd()}/setup")
        self.pdb_dir = pathlib.Path(f'{self.setup_rundir}/analysis')
        self.Setup = Setup(f'{data_folder}/setup.toml')

    def test_initialize(self):
        observed_run_params = self.Setup.initialize()
        expected_run_params = {'folder': f'{self.setup_rundir}',
                               'mol_a': f'{self.setup_rundir}/input/molA.pdb',
                               'mol_b': f'{self.setup_rundir}/input/molB.pdb',
                               'np': 1,
                               'restraints_a': [39, 40, 41],
                               'restraints_b': [4, 5, 6]}
        self.assertEqual(observed_run_params, expected_run_params)

    def test_clean(self):
        self.pdb_dir.mkdir(parents=True)
        expected_directory_contents = []
        for i in range(10):
            shutil.copy(f'{data_folder}/molA.pdb', f'{self.pdb_dir}/{i}.pdb')
            with open(f'{self.pdb_dir}/{i}.contacts', 'w') as con_fh:
                con_fh.write('11111111')
            con_fh.close()
            expected_directory_contents.append(f'{self.pdb_dir}/{i}.pdb.gz')

        self.Setup.clean()
        observed_directory_contents = [x for x in glob.glob(f'{self.pdb_dir}/*')]
        observed_directory_contents.sort()

        self.assertListEqual(observed_directory_contents, expected_directory_contents)

    def tearDown(self):
        shutil.rmtree(pathlib.Path(f"{pathlib.Path.cwd()}/setup"))
