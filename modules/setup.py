import os
import shutil
import toml
import logging
import glob
import sys
import configparser
import pathlib
import subprocess
import tempfile
import shlex
import multiprocessing
import mgzip
from utils.files import get_full_path
from utils.functions import du, check_if_py3
from modules.error import (DependencyNotDefinedError, DependencyNotFoundError,
                           SectionNotDefinedError)
ga_log = logging.getLogger('ga_log')


class Setup:
    def __init__(self, toml_file):
        """Initialize setup class."""
        self.input_params = toml.load(toml_file)
        self.data = {}
        self.etc_folder = get_full_path('etc')
        self.nproc = self.input_params['main']['number_of_processors']

    def initialize(self):
        """Load the parameters and create the folder structure."""
        ga_params = toml.load(f"{self.etc_folder}/"
                              "genetic_algorithm_params.toml")
        run_params = {}
        identifier_folder = self.input_params['main']['identifier']
        run_path = f'{os.getcwd()}/{identifier_folder}'
        ga_log.info('Initializing')

        ga_log.info('Checking for dependencies')
        try:
            self.validate_third_party()
        except Exception as e:
            ga_log.error(e)
            sys.exit()

        ga_log.debug(f'Run path: {run_path}')
        ga_log.debug(f'Run folder: {identifier_folder}')

        if os.path.isdir(identifier_folder):
            ga_log.warning(f'Your run folder {identifier_folder}'
                           ' will be deleted!')
            shutil.rmtree(identifier_folder)

        os.mkdir(identifier_folder)

        mol_a = self.input_params['molecules']['A']
        mol_a_name = mol_a.split('/')[-1]
        mol_b = self.input_params['molecules']['B']
        mol_b_name = mol_b.split('/')[-1]
        if 'native' in self.input_params['molecules']:
            native = self.input_params['molecules']['native']
        else:
            native = ''
        input_folder = f'{identifier_folder}/input'

        ga_log.info('Copying input molecules to run folder')
        if not os.path.isdir(input_folder):
            ga_log.debug(f'Creating input folder: {input_folder}')
            os.mkdir(input_folder)
            ga_log.debug(f'Copying {mol_a}')
            shutil.copy(mol_a, input_folder)
            ga_log.debug(f'Copying {mol_b}')
            shutil.copy(mol_b, input_folder)
            if native:
                ga_log.debug(f'Copying {native}')
                shutil.copy(native, input_folder)
                izone = native.replace('.pdb', '.izone')
                izone = pathlib.Path(izone)
                if izone.exists():
                    shutil.copy(izone, input_folder)

        analysis_folder = f'{identifier_folder}/analysis'
        if not os.path.isdir(analysis_folder):
            ga_log.debug(f'Creating analysis folder {analysis_folder}')
            os.mkdir(analysis_folder)

        structures_folder = f'{identifier_folder}/structures'
        if not os.path.isdir(structures_folder):
            ga_log.debug(f'Creating structures folder {structures_folder}')
            os.mkdir(structures_folder)

        run_params['folder'] = run_path
        run_params['mol_a'] = f'{run_path}/input/{mol_a_name}'
        run_params['mol_b'] = f'{run_path}/input/{mol_b_name}'
        if native:
            run_params['native'] = f'{run_path}/input/{native}'
        run_params['restraints_a'] = self.input_params['restraints']['A']
        run_params['restraints_b'] = self.input_params['restraints']['B']
        run_params['np'] = self.input_params['main']['number_of_processors']

        return run_params, ga_params

    def clean(self):
        """Clean the run directory."""
        ga_log.info('Cleaning the simulation directory')
        identifier_folder = self.input_params['main']['identifier']

        run_path = f'{os.getcwd()}/{identifier_folder}'
        structure_path = f'{os.getcwd()}/{identifier_folder}/structures'
        analysis_path = f'{os.getcwd()}/{identifier_folder}/analysis'
        fcc_matrix_f = f'{analysis_path}/fcc.matrix'

        pool = multiprocessing.Pool(processes=self.nproc)

        size = du(run_path)
        ga_log.info(f'Current size: {size}')

        ga_log.info('Compressing PDB structures')
        pdb_list = glob.glob(f'{structure_path}/*pdb')
        pool.map_async(self.compress, pdb_list)

        ga_log.info('Deleting .contacts files')
        contact_list = glob.glob(f'{structure_path}/*contacts')
        pool.map_async(os.remove, contact_list)

        pool.close()
        pool.join()

        if os.path.isfile(fcc_matrix_f):
            ga_log.info('Compressing fcc.matrix')
            self.compress(fcc_matrix_f, np=self.nproc)

        size = du(run_path)
        ga_log.info(f'Cleaning done - current size: {size}')

    @staticmethod
    def compress(file_path, np=1):
        """Compress a file in .gz."""
        # Note: this will open the whole file in memory
        #  this might not be the best idea
        with open(file_path, 'r') as fr:
            file_string = ''.join(fr.readlines())
        fr.close()

        with mgzip.open(f'{file_path}.gz', "wt", thread=np) as fw:
            fw.write(file_string)
        fw.close()

        os.remove(file_path)

    def validate_third_party(self):
        """Check if the third-party dependencies are ok in gdock.ini."""
        ga_ini = configparser.ConfigParser(os.environ)
        ga_ini.read(os.path.join(self.etc_folder, 'gdock.ini'),
                    encoding='utf-8')

        try:
            ga_ini.get('third_party', '')
        except configparser.NoSectionError:
            raise SectionNotDefinedError('third_party')
        except configparser.NoOptionError:
            pass

        # Check FCC
        try:
            fcc_path = pathlib.Path(ga_ini.get('third_party', 'fcc_path'))
        except configparser.NoOptionError:
            raise DependencyNotDefinedError('fcc_path')

        if not fcc_path.exists():
            raise DependencyNotFoundError(f'{fcc_path} not found')

        fcc_script = pathlib.Path(f'{fcc_path}/scripts/make_contacts.py')
        if fcc_script.exists():
            fcc_check = check_if_py3(str(fcc_script))
            if not fcc_check:
                raise Exception('FCC branch is not Python3 compatible,'
                                ' check install instructions')
        else:
            raise Exception(f'{fcc_script} not found')

        # Check pdbtools
        try:
            pdbtools_path = pathlib.Path(ga_ini.get('third_party',
                                                    'pdbtools_path'))
        except configparser.NoOptionError:
            raise DependencyNotDefinedError('pdbtools_path')

        if not pdbtools_path.exists():
            raise DependencyNotFoundError(f'{pdbtools_path} not found')

        pdbtools_script = pathlib.Path(f'{pdbtools_path}/pdbtools/pdb_tidy.py')
        if pdbtools_script.exists():
            pdbtools_check = check_if_py3(str(pdbtools_script))
            if not pdbtools_check:
                raise Exception('pdb-tools branch is not Python3 compatible,'
                                ' check install instructions')
        else:
            raise Exception(f'{pdbtools_script} not found')

        # check haddock-tools
        try:
            haddocktools_path = pathlib.Path(ga_ini.get('third_party',
                                                        'haddocktools_path'))
        except configparser.NoOptionError:
            raise DependencyNotDefinedError('haddocktools_path')

        if not haddocktools_path.exists():
            raise DependencyNotFoundError(f'{haddocktools_path} not found')

        haddocktools_script = pathlib.Path((f'{haddocktools_path}/'
                                            'contact-chainID'))
        if haddocktools_script.exists():
            # check if executable
            proc = subprocess.run(str(haddocktools_script),
                                  stderr=subprocess.PIPE,
                                  stdout=subprocess.PIPE)

            err = proc.stderr.decode('utf-8')
            out = proc.stdout.decode('utf-8')
            if 'Too few arguments' not in err:
                raise Exception(f'{haddocktools_script} execution failed', err)

        else:
            raise Exception(f'{haddocktools_script} not found')

        # Check dcomplex
        try:
            dcomplex_exe = pathlib.Path(ga_ini.get('third_party',
                                                   'dcomplex_exe'))
        except configparser.NoOptionError:
            raise DependencyNotDefinedError('dcomplex_exe')

        if dcomplex_exe.exists():
            # check if executable
            proc = subprocess.run(str(dcomplex_exe),
                                  stderr=subprocess.PIPE,
                                  stdout=subprocess.PIPE)

            err = proc.stderr.decode('utf-8')
            out = proc.stdout.decode('utf-8')
            if 'usage' not in out:
                raise Exception(f'{dcomplex_exe} execution failed', err)
        else:
            raise DependencyNotFoundError(dcomplex_exe)

        # Check PROFIT
        try:
            profit_exe = pathlib.Path(ga_ini.get('third_party',
                                                 'profit_exe'))
        except configparser.NoOptionError:
            raise DependencyNotDefinedError('profit_exe')

        if profit_exe.exists():
            # check if executable
            dummy_script = tempfile.NamedTemporaryFile(delete=False,
                                                       suffix='.txt')
            cmd = f'{profit_exe} -f {dummy_script.name}'
            out = subprocess.check_output(shlex.split(cmd), shell=False)
            out = out.decode('utf-8')
            os.unlink(dummy_script.name)

            if 'Finished script' not in out:
                raise Exception(f'{profit_exe} execution failed')
        else:
            raise DependencyNotFoundError(profit_exe)

        return True
