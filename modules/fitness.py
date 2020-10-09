from subprocess import Popen, PIPE
import os
import configparser
import numpy as np
from utils.files import get_full_path
from utils.functions import get_coords
import logging
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
    proc = Popen(cmd, shell=True, stdout=PIPE)
    result = proc.stdout.read().decode('utf-8')
    irmsd = float(result.split('\n')[-6].split()[-1])
    return irmsd
