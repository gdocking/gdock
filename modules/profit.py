import configparser
import os
import subprocess
import shlex
import multiprocessing
import pathlib
import numpy as np
import logging
from tempfile import NamedTemporaryFile
from utils.files import get_full_path

ga_log = logging.getLogger('ga_log')
etc_folder = get_full_path('etc')
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, 'gdock.ini'), encoding='utf-8')
profit_exe = ini.get('third_party', 'profit_exe')


class Profit:
    """Wrapper for PROFIT."""

    def __init__(self, ref, mobi, nproc):
        self.reference = ref  # str
        self.mobi = mobi  # list
        self.exec = profit_exe
        self.nproc = int(nproc)
        self.izone = ''
        izone_f = pathlib.Path(ref.replace('.pdb', '.izone'))
        if not izone_f.exists():
            raise Exception(f'{izone_f} not found')
        else:
            with open(izone_f, 'r') as fh:
                for line in fh.readlines():
                    # the izones come with \n from BM5,
                    #  so split them here explicitly
                    line = line.split('\n')[0]
                    self.izone += line + os.linesep

    def calc_irmsd(self):
        """Calculate the irmsd using PROFIT."""
        irmsd_dic = {}
        pool = multiprocessing.Pool(processes=self.nproc)

        results = []
        for structure in self.mobi:
            # Q: how to clean the temporary scripts?
            script = self._write_script(structure)
            results.append((structure,
                            pool.apply_async(self.execute,
                                             args=(self.exec, script.name)))
                           )

        for p in results:
            pdb, irmsd = p[0], p[1].get()
            pdb_identifier = pathlib.Path(pdb).name[:-4]
            irmsd_dic[pdb_identifier] = irmsd

        pool.close()
        pool.join()

        # do some analysis
        irmsd_list = list(irmsd_dic.values())
        irmsd_quantiles = np.quantile(irmsd_list, [0.0, 0.25, 0.5, 0.75, 1.0])
        irmsd_mean = irmsd_quantiles[2]
        irmsd_min = irmsd_quantiles[0]
        irmsd_max = irmsd_quantiles[4]
        irmsd_sd = np.std(irmsd_list)
        ga_log.info(f'min: {irmsd_min:.2f} Å max: {irmsd_max:.2f} Å mean: '
                    f'{irmsd_mean:.2f} Å sd: {irmsd_sd:.2f} Å')

        return irmsd_dic

    def _write_script(self, mobi):
        """Write a script to a temporary file for PROFIT."""
        script_str = f'ref {self.reference}' + os.linesep
        script_str += f'mobi {mobi}' + os.linesep
        script_str += 'atoms C,CA,N,O' + os.linesep
        script_str += self.izone
        script_str += 'fit' + os.linesep
        script_str += 'quit'

        script_f = NamedTemporaryFile(delete=False, suffix='.txt')
        script_f.write(str.encode(script_str))
        script_f.close()

        return script_f

    @staticmethod
    def execute(profit_exe, script_f):
        """Execute PROFIT."""
        cmd = f'{profit_exe} -f {script_f}'
        try:
            out = subprocess.check_output(shlex.split(cmd), shell=False)
            result = out.decode('utf-8')
            irmsd = float(result.split()[-1])
        except Exception as e:
            ga_log.warning(e)
            irmsd = float('nan')

        return irmsd
