import subprocess
import numpy as np

dcomplex_exe = '/Users/rodrigo/repos/gadock/src/dcomplex_single_file/dcomplex'
dockq_exe = '/Users/rodrigo/repos/DockQ/DockQ.py'

def calc_irmsd(pdb_f):
    native = '/Users/rodrigo/repos/gadock/target-unbound.pdb'
    cmd = f'{dockq_exe} {pdb_f} {native}'
    proc = subprocess.Popen(cmd.split(), shell=False, stdout=subprocess.PIPE)
    result = proc.stdout.read().decode('utf-8')
    irmsd = float(result.split('\n')[-6].split()[-1])
    return irmsd

def dcomplex(pdb_f):
    cmd = f'{dcomplex_exe} {pdb_f} A B'
    proc = subprocess.Popen(cmd.split(), shell=False, stdout=subprocess.PIPE)
    energ = float(proc.stdout.read().decode('utf-8').split('\n')[-2].split()[1])
    return energ

def calc_clash(input_pdb):
    d = {}
    with open(input_pdb, 'r') as fh:
        for l in fh.readlines():
            if l.startswith('ATOM'):
                chain = l[21]
                if chain not in d:
                    d[chain] = []
                x = float(l[31:38])
                y = float(l[39:46])
                z = float(l[47:54])
                d[chain].append((x,y,z))
    #
    distances_list = []
    for chain_x in d:
        for coord_a in d[chain_x]:
            xa, ya, za  = coord_a
            for chain_y in d:
                if chain_x != chain_y:
                    for coord_b in d[chain_y]:
                        xb, yb, zb = coord_b
                        dist = np.sqrt((xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2)
                        distances_list.append(dist)
    clash_score = len([d for d in distances_list if d <= 4.0])
    return clash_score