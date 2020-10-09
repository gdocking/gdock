import os
import shutil
import toml
import logging
ga_log = logging.getLogger('ga_log')


class Setup:
    def __init__(self, toml_file):
        """Initialize setup class.


        :param toml_file:
        """
        self.input_params = toml.load(toml_file)
        self.data = {}

    def initialize(self):
        """Load the parameters and create the folder structure."""
        run_params = {}
        identifier_folder = self.input_params['main']['identifier']
        run_path = f'{os.getcwd()}/{identifier_folder}'
        ga_log.info('Initializing')
        ga_log.debug(f'Run path: {run_path}')
        ga_log.debug(f'Run folder: {identifier_folder}')
        if not os.path.isdir(identifier_folder):
            os.mkdir(identifier_folder)

        mol_a = self.input_params['molecules']['A'].split('/')[-1]
        mol_b = self.input_params['molecules']['B'].split('/')[-1]
        input_folder = f'{identifier_folder}/input'

        ga_log.info('Copying input molecules to run folder')
        if not os.path.isdir(input_folder):
            ga_log.debug(f'Creating input folder: {input_folder}')
            os.mkdir(input_folder)
            ga_log.debug(f'Copying {mol_a}')
            shutil.copy(mol_a, input_folder)
            ga_log.debug(f'Copying {mol_b}')
            shutil.copy(mol_b, input_folder)

        begin_folder = f'{identifier_folder}/begin'
        if not os.path.isdir(begin_folder):
            ga_log.debug(f'Creating begin folder {begin_folder}')
            os.mkdir(begin_folder)

        gen_folder = f'{identifier_folder}/gen'
        if not os.path.isdir(gen_folder):
            ga_log.debug(f'Creating gen folder {gen_folder}')
            os.mkdir(gen_folder)

        run_params['folder'] = run_path
        run_params['mol_a'] = f'{run_path}/input/{mol_a}'
        run_params['mol_b'] = f'{run_path}/input/{mol_b}'
        run_params['restraints_a'] = self.input_params['restraints']['A']
        run_params['restraints_b'] = self.input_params['restraints']['B']
        run_params['np'] = self.input_params['main']['number_of_processors']

        return run_params
