# Deamon to run the benchmark!
import argparse
import glob
import logging
import os
import shlex
import subprocess  # nosec
import sys
from pathlib import Path

bm_log = logging.getLogger("bm_log")
bm_log.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
formatter = logging.Formatter(
    " %(asctime)s %(module)s:%(lineno)d" " %(levelname)s - %(message)s"
)
ch.setFormatter(formatter)
bm_log.addHandler(ch)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "gdockbm_path",
        help=("Full path of the root directory" " of the prepared folders"),
    )
    args = parser.parse_args()

    benchmark_path = Path(args.gdockbm_path).resolve()
    if not benchmark_path.exists():
        bm_log.error(f"Path {benchmark_path} does not exist.")
        sys.exit()

    # FIXME: There is probably a better way of doing this
    dirname, filename = os.path.split(os.path.abspath(__file__))
    gdock_exe = f"{dirname}/../gdock.py"
    python_exe = sys.executable

    benchmark = [f for f in glob.glob(f"{benchmark_path}/*") if "." not in f]
    benchmark.sort()

    total = len(benchmark)

    for i, folder in enumerate(benchmark):
        target_name = Path(folder).name
        output_file = f"{folder}/run/analysis/gdock.dat"
        if not os.path.isfile(output_file):

            bm_log.info(f"{target_name} Running - Remaining {total-i}")

            os.chdir(folder)
            bm_log.debug(f"chdir {folder}")

            cmd = f"{python_exe} {gdock_exe} run.toml"
            bm_log.debug(f"cmd is: {cmd}")
            result = subprocess.run(  # nosec
                shlex.split(cmd), capture_output=True, shell=False
            )

            # Make this explicit so eventually we can capture
            #  stuff from stderr/out
            stdout = result.stdout.decode("utf-8")
            stderr = result.stderr.decode("utf-8")

            with open(f"{folder}/gdock.out", "w") as out_fh:
                out_fh.write(stderr)
            out_fh.close()

            if not os.path.isfile(output_file):
                bm_log.error(f"Something wrong with {target_name}")
                with open(f"{folder}/FAIL", "w") as fail_fh:
                    fail_fh.write("I failed :(")
                fail_fh.close()

        else:
            bm_log.info(f"{target_name} DONE - Remaining {total-i}")
