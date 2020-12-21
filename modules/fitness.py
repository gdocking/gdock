import os
import logging
import configparser
import subprocess  # nosec
import shlex
from utils.files import get_full_path
ga_log = logging.getLogger('ga_log')

etc_folder = get_full_path('etc')
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, 'gdock.ini'), encoding='utf-8')
dcomplex_exe = ini.get('third_party', 'dcomplex_exe')


def run_dcomplex(pdb_f):
    """Use DCOMPLEX to calculate the PDBs energy."""
    cmd = (f'{dcomplex_exe} {pdb_f} A B')
    ga_log.debug(f'cmd is: {cmd}')
    out = subprocess.check_output(shlex.split(cmd), shell=False)  # nosec
    result = out.decode('utf-8').split('\n')
    energy = float(result[-2].split()[1])
    return energy
