# Deamon to run the benchmark!
import glob
import os
import sys
from pathlib import Path
import logging
import argparse
import subprocess  # nosec
import shlex

bm_log = logging.getLogger('bm_log')
bm_log.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
formatter = logging.Formatter(' %(asctime)s %(module)s:%(lineno)d %(levelname)s - %(message)s')
ch.setFormatter(formatter)
bm_log.addHandler(ch)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("gdockbm_path", help="Location where the prepared folders will be")
    args = parser.parse_args()

    # FIXME: There is probably a better way of doing this
    dirname, filename = os.path.split(os.path.abspath(__file__))
    gdock_exe = f'{dirname}/../gdock.py'
    python_exe = sys.executable

    benchmark = glob.glob(f'{args.gdockbm_path}/*')
    benchmark.sort()

    total = len(benchmark)

    for i, folder in enumerate(benchmark):
        target_name = Path(folder).name
        bm_log.info(f'{target_name} Running - Remaining {total-i}')

        os.chdir(folder)
        bm_log.debug(f'chdir {folder}')

        cmd = (f'{python_exe} {gdock_exe} run.toml')
        bm_log.debug(f'cmd is: {cmd}')
        result = subprocess.run(shlex.split(cmd), capture_output=True, shell=False)  # nosec

        # Make this explicit so eventually we can capture stuff from stderr/out
        stdout = result.stdout.decode('utf-8')
        stderr = result.stderr.decode('utf-8')

        with open(f'{folder}/gdock.out', 'w') as out_fh:
            out_fh.write(stderr)

        if result.returncode != 0:
            bm_log.warning(f'Something went wrong with {target_name}')
        elif result.returncode == 0:
            bm_log.info(f'{target_name} Complete')

        # FIXME: WIP
        if i == 5:
            sys.exit()
