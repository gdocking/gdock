# this script will read the finished benchmark and prepare a file for analysis
import logging
import pathlib
import argparse
import os

bm_log = logging.getLogger('bm_log')
bm_log.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
formatter = logging.Formatter(' %(asctime)s %(module)s:%(lineno)d'
                              ' %(levelname)s - %(message)s')
ch.setFormatter(formatter)
bm_log.addHandler(ch)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("bm_path", help=("Full path of the root directory"
                                         " of the prepared folders"))
    parser.add_argument("output_f", help=("Name of the output benchmark"
                                          " data file"))
    args = parser.parse_args()

    benchmark_path = pathlib.Path(args.bm_path)
    benchmark_f = pathlib.Path(args.output_f)

    bm_log.info('Retrieving benchmark data')
    bm_log.info(f'Benchmark path is {benchmark_path}')
    bm_log.info(f'Output file is {benchmark_f}')

    # get the header
    for f in benchmark_path.glob('*/run/analysis/gdock.dat'):
        header = open(f).readlines()[0].split()
        header.append('target_name')
        bm_log.info(f'Header is {",".join(header)}')
        break

    with open(benchmark_f, 'w') as data_fh:
        data_fh.write('\t'.join(header) + os.linesep)
        for data_path in benchmark_path.glob('*/run/analysis/gdock.dat'):
            target_name = str(data_path).split('/')[-4]
            with open(data_path) as in_fh:
                for line in in_fh.readlines():
                    data = line.split()
                    if 'gen' not in data:
                        data.append(target_name)
                        new_line = '\t'.join(data) + os.linesep
                        data_fh.write(new_line)
            in_fh.close()
    data_fh.close()
