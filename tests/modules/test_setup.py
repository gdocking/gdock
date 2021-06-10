import unittest
import pathlib
import shutil
import glob
import os
import tempfile
import copy
from modules.setup import Setup
from modules.error import (DependencyNotDefinedError, DependencyNotFoundError,
                           SectionNotDefinedError)
from utils.files import get_full_path

data_folder = get_full_path('tests', 'test_data')


class TestSetup(unittest.TestCase):
    def setUp(self):

        toml_string = "[main]" + os.linesep
        toml_string += "identifier = 'setup'" + os.linesep
        toml_string += "number_of_processors = 1" + os.linesep
        toml_string += "" + os.linesep
        toml_string += "[restraints]" + os.linesep
        toml_string += "A = [39,40,41]" + os.linesep
        toml_string += "B = [4,5,6]" + os.linesep
        toml_string += "" + os.linesep
        toml_string += "[molecules]" + os.linesep
        toml_string += f"A = '{data_folder}/molA.pdb'" + os.linesep
        toml_string += f"B = '{data_folder}/molB.pdb'" + os.linesep

        with open(f'{data_folder}/setup.toml', 'w') as setup_f:
            setup_f.write(toml_string)
        self.setup_rundir = pathlib.Path(f"{pathlib.Path.cwd()}/setup")
        self.pdb_dir = pathlib.Path(f'{self.setup_rundir}/structures')
        self.analysis_dir = pathlib.Path(f'{self.setup_rundir}/analysis')
        self.Setup = Setup(f'{data_folder}/setup.toml')

    def test_initialize(self):
        observed_run_params, _ = self.Setup.initialize()
        expected_run_params = {'folder': f'{self.setup_rundir}',
                               'mol_a': f'{self.setup_rundir}/input/molA.pdb',
                               'mol_b': f'{self.setup_rundir}/input/molB.pdb',
                               'np': 1,
                               'restraints_a': [39, 40, 41],
                               'restraints_b': [4, 5, 6]}
        self.assertEqual(observed_run_params, expected_run_params)

    def test_clean(self):
        """Test final cleaning."""
        # create dummy folder
        self.pdb_dir.mkdir(exist_ok=True, parents=True)
        # create dummy pdbs
        for i in range(10):
            shutil.copy(f'{data_folder}/molA.pdb', f'{self.pdb_dir}/{i}.pdb')

        # create dummy contact files
        for i in range(10):
            with open(f'{self.pdb_dir}/{i}.contacts', 'w') as con_fh:
                con_fh.write('11111111')
            con_fh.close()

        # create dummy fcc.matrix
        self.analysis_dir.mkdir(exist_ok=True, parents=True)
        with open(f'{self.analysis_dir}/fcc.matrix', 'w') as out:
            out.seek((1024 * 1024) - 1)
            out.write('\0')
        out.close()

        self.Setup.clean()

        expected_pdb_f = [f'{self.pdb_dir}/{i}.pdb.gz' for i in range(10)]
        observed_pdb_f = [x for x in glob.glob(f'{self.pdb_dir}/*pdb*')]

        observed_pdb_f.sort()
        expected_pdb_f.sort()
        self.assertListEqual(observed_pdb_f, expected_pdb_f)

        expected_con_f = []
        observed_con_f = [x for x in glob.glob(f'{self.pdb_dir}/*con*')]
        self.assertListEqual(observed_con_f, expected_con_f)

        expected_m_f = pathlib.Path(f'{self.analysis_dir}/fcc.matrix.gz')
        self.assertTrue(expected_m_f.exists())

        expected_m_f_size = 1072
        observed_m_f_size = expected_m_f.stat().st_size

        self.assertEqual(observed_m_f_size, expected_m_f_size)

    def test_compress(self):
        """Test file compression."""

        dummy_f = tempfile.NamedTemporaryFile(delete=False)
        dummy_f.seek((1024 * 1024) - 1)
        dummy_f.write(b'\0')
        dummy_f.close()

        self.Setup.compress(dummy_f.name, np=1)

        compressed_file = pathlib.Path(dummy_f.name + '.gz')

        self.assertTrue(compressed_file.exists())

        expected_compressed_fsize = 1073
        observed_compressed_fsize = compressed_file.stat().st_size

        self.assertEqual(observed_compressed_fsize, expected_compressed_fsize)

        compressed_file.unlink()

    def test_validate_third_party(self):
        """Test third-party validation."""
        valid_ini = copy.deepcopy(self.Setup.ga_ini)

        self.Setup.ga_ini.remove_section('third_party')
        self.assertRaises(SectionNotDefinedError,
                          self.Setup.validate_third_party)
        self.Setup.ga_ini = copy.deepcopy(valid_ini)

        self.Setup.ga_ini.set('third_party', 'fcc_path', 'not_valid')
        self.assertRaises(DependencyNotFoundError,
                          self.Setup.validate_third_party)
        self.Setup.ga_ini = copy.deepcopy(valid_ini)

        self.Setup.ga_ini.remove_option('third_party', 'pdbtools_path')
        self.assertRaises(DependencyNotDefinedError,
                          self.Setup.validate_third_party)
        self.Setup.ga_ini = valid_ini

    def tearDown(self):
        if self.setup_rundir.exists():
            shutil.rmtree(self.setup_rundir)
