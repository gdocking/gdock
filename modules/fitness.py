import os
import logging
import configparser
import subprocess
import shlex
from utils.files import get_full_path
ga_log = logging.getLogger('ga_log')

etc_folder = get_full_path('etc')
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, 'gadock.ini'), encoding='utf-8')
native = ini.get('reference', 'native_pdb')

dockq_exe = ini.get('third_party', 'dockq_exe')


def calc_irmsd(pdb_f):
    """Calculate the interface root mean square deviation."""
    cmd = f'{dockq_exe} {pdb_f} {native}'
    ga_log.debug(f'cmd is: {cmd}')
    out = subprocess.check_output(shlex.split(cmd), shell=False)
    result = out.decode('utf-8')
    irmsd = float(result.split('\n')[-6].split()[-1])
    return irmsd
