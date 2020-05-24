#PBS -N GADock
#PBS -q medium
#PBS -l nodes=1:ppn=1
#PBS -S /bin/tcsh

cd /home/rodrigo/repos/gadock
setenv PYTHONPATH ${PYTHONPATH}:`pwd`
conda activate gadock
python bin/gadock.py target-unbound.pdb 10 25 1 > gadock.out
