#PBS -N GADock
#PBS -q medium
#PBS -l nodes=1:ppn=48
#PBS -S /bin/tcsh

cd /home/rodrigo/repos/gadock
setenv PYTHONPATH ${PYTHONPATH}:`pwd`
conda activate gadock
python bin/gadock.py target-unbound.pdb --np 48 > gadock.out
