"""Setup Module."""
import os
import shutil
import logging
import toml
ga_log = logging.getLogger('ga_log')


class Setup:
    """Setup Class."""

    def __init__(self, toml_file):
        """Initialize setup class."""
        self.input_params = toml.load(toml_file)
        self.data = {}

    def initialize(self):
        """Load the parameters and create the folder structure."""
        run_params = {}
        identifier_folder = self.input_params['main']['identifier']
        run_path = f'{os.getcwd()}/{identifier_folder}'
        ga_log.info('Initializing')
        ga_log.debug('Run path: %s', run_path)
        ga_log.debug('Run folder: %s', identifier_folder)
        if not os.path.isdir(identifier_folder):
            os.mkdir(identifier_folder)

        mol_a = self.input_params['molecules']['A'].split('/')[-1]
        mol_b = self.input_params['molecules']['B'].split('/')[-1]
        input_folder = f'{identifier_folder}/input'

        ga_log.info('Copying input molecules to run folder')
        if not os.path.isdir(input_folder):
            ga_log.debug('Creating input folder: %s', input_folder)
            os.mkdir(input_folder)
            ga_log.debug('Copying %s', mol_a)
            shutil.copy(mol_a, input_folder)
            ga_log.debug('Copying %s', mol_b)
            shutil.copy(mol_b, input_folder)

        begin_folder = f'{identifier_folder}/begin'
        if not os.path.isdir(begin_folder):
            ga_log.debug('Creating begin folder %s', begin_folder)
            os.mkdir(begin_folder)

        gen_folder = f'{identifier_folder}/gen'
        if not os.path.isdir(gen_folder):
            ga_log.debug('Creating gen folder %s', gen_folder)
            os.mkdir(gen_folder)

        analysis_folder = f'{identifier_folder}/analysis'
        if not os.path.isdir(analysis_folder):
            ga_log.debug('Creating analysis folder %s', gen_folder)
            os.mkdir(analysis_folder)

        run_params['folder'] = run_path
        run_params['mol_a'] = f'{run_path}/input/{mol_a}'
        run_params['mol_b'] = f'{run_path}/input/{mol_b}'
        run_params['restraints_a'] = self.input_params['restraints']['A']
        run_params['restraints_b'] = self.input_params['restraints']['B']
        run_params['np'] = self.input_params['main']['number_of_processors']

        return run_params
