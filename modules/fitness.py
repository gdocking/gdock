import os
import logging
import configparser
import subprocess  # nosec
import shlex
import scipy
from utils.files import get_full_path
ga_log = logging.getLogger('ga_log')

etc_folder = get_full_path('etc')
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, 'gdock.ini'), encoding='utf-8')
dcomplex_exe = ini.get('third_party', 'dcomplex_exe')

# ============================================= #
# Keep the fitness functions self-contained!    #
# ============================================= #


def run_dcomplex(pdb_f):
    """Use DCOMPLEX to calculate the PDB's energy."""
    cmd = (f'{dcomplex_exe} {pdb_f} A B')
    ga_log.debug(f'cmd is: {cmd}')
    out = subprocess.check_output(shlex.split(cmd), shell=False)  # nosec
    result = out.decode('utf-8').split('\n')
    energy = float(result[-2].split()[1])
    return energy


# TODO: speed this up!
def calc_satisfaction(pdb_f, restraints_a, restraints_b, cutoff=4.9):
    """Calculate the restraints satisfaction ratio."""
    coord_dic = {}
    with open(pdb_f) as fh:
        for line in fh:
            if line.startswith('ATOM'):
                chain = line[21]
                if chain not in coord_dic:
                    coord_dic[chain] = {}
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                resnum = int(line[22:26])

                if chain == 'A' and resnum in restraints_a:
                    if resnum not in coord_dic[chain]:
                        coord_dic[chain][resnum] = []
                    coord_dic[chain][resnum].append((x, y, z))

                if chain == 'B' and resnum in restraints_b:
                    if resnum not in coord_dic[chain]:
                        coord_dic[chain][resnum] = []
                    coord_dic[chain][resnum].append((x, y, z))

    coord_a = coord_dic['A']
    coord_b = coord_dic['B']

    contacts = {'A': [], 'B': []}
    for res_a in coord_a:
        for atom_coords_a in coord_a[res_a]:
            x_a, y_a, z_a = atom_coords_a
            for res_b in coord_b:
                for atom_coords_b in coord_b[res_b]:
                    x_b, y_b, z_b = atom_coords_b

                    c_i = x_a, y_a, z_a
                    c_j = x_b, y_b, z_b

                    distance = scipy.spatial.distance.euclidean(c_i, c_j)

                    if distance <= cutoff:
                        contacts['A'].append(res_a)
                        contacts['B'].append(res_b)

    # check how many are satisfied
    satisfied_a = len(set(contacts['A']).intersection(restraints_a))
    satisfied_b = len(set(contacts['B']).intersection(restraints_b))

    total_restraints = len(restraints_a) + len(restraints_b)

    satisfied_ratio = (satisfied_a + satisfied_b) / total_restraints

    return satisfied_ratio
