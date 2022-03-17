#!/usr/bin/bash
#SBATCH --partition=verylong
#SBATCH --ntasks-per-node=48
#SBATCH --nodes=1
#SBATCH --job-name="gdock_bm5"
#SBATCH --output=bm.out

BM_PATH=

source $HOME/software/miniconda3/etc/profile.d/conda.sh
conda activate gdock
cd BM_PATH
python $HOME/repos/gdock/tools/run-benchmark.py BM_PATH >& bm.out