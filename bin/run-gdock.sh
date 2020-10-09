#!/usr/bin/env bash
#PBS -N gdock
#PBS -q medium
#PBS -l nodes=1:ppn=48
#PBS -S /bin/tcsh

cd /home/rodrigo/repos/gdock
setenv "PYTHONPATH ${PYTHONPATH}:$(pwd)"
conda activate gadock
python bin/gdock.py --np 48 > gdock.out
