import os
import logging
import configparser
import subprocess  # nosec
import shlex
import tempfile
from utils.files import get_full_path
from pathlib import Path
ga_log = logging.getLogger('ga_log')

etc_folder = get_full_path('etc')
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, 'gdock.ini'), encoding='utf-8')
foldx_exe = ini.get('third_party', 'foldx_exe')


def run_foldx(pdb_f):
    """Use FoldX to calculate the PDBs energy."""
    # in FoldX - More negative energies indicate better binding. Positive energies indicate no binding.
    # foldx needs both the pdbfilename and where it is
    loc = Path(pdb_f)
    pdb_location = str(loc.parent)
    pdb_filename = loc.name
    # if we dont pass a temporary directory it will write some files on the running dir
    temp_dir = tempfile.TemporaryDirectory()
    cmd = (f'{foldx_exe} -c AnalyseComplex --pdb {pdb_filename} --noHeader 1 --pdb-dir={pdb_location} '
           f'--output-dir={temp_dir.name}')
    ga_log.debug(f'cmd is: {cmd}')
    try:
        out = subprocess.check_output(shlex.split(cmd), shell=False)  # nosec
        result = out.decode('utf-8').split('\n')
        # FIXME: put some fancy regex here
        energy = float(result[-8].split()[-1])
    except subprocess.CalledProcessError as e:
        ga_log.error(f'Foldx failed with {e}')
        energy = float('nan')
    temp_dir.cleanup()
    return energy
