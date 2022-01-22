import configparser
import logging
import os
import shlex
import subprocess  # nosec

# from utils.functions import timer
from gdock.modules.files import get_full_path

ga_log = logging.getLogger("ga_log")

etc_folder = get_full_path("etc")
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, "gdock.ini"), encoding="utf-8")
haddocktools_path = ini.get("third_party", "haddocktools_path")

# ============================================= #
# Keep the fitness functions self-contained!    #
# ============================================= #


def calc_satisfaction(pdb_f, restraints_a, restraints_b, cutoff=4.9):
    """Calculate the restraints satisfaction ratio."""
    # this is 4x faster!
    cmd = f"{haddocktools_path}/contact-chainID {pdb_f} {cutoff}"
    out = subprocess.check_output(shlex.split(cmd))  # nosec
    contacts = {"A": [], "B": []}
    for line in out.decode("utf-8").split(os.linesep):
        data = line.split()
        if data:
            res_a, chain_a, _, res_b, chain_b, _, _ = data
            res_a = int(res_a)
            res_b = int(res_b)
            if chain_a == "A" and res_a in restraints_a:
                contacts["A"].append(res_a)
            if chain_b == "B" and res_b in restraints_b:
                contacts["A"].append(res_a)

    # check how many are satisfied
    satisfied_a = len(set(contacts["A"]).intersection(restraints_a))
    satisfied_b = len(set(contacts["B"]).intersection(restraints_b))

    total_restraints = len(restraints_a) + len(restraints_b)

    satisfied_ratio = (satisfied_a + satisfied_b) / total_restraints

    return satisfied_ratio


# TODO: Do a low-level implementation of this
# def calc_satisfaction(pdb_f, restraints_a, restraints_b, cutoff=4.9):
#     """Calculate the restraints satisfaction ratio."""
#     coord_dic = {}
#     with open(pdb_f) as fh:
#         for line in fh:
#             if line.startswith('ATOM'):
#                 chain = line[21]
#                 if chain not in coord_dic:
#                     coord_dic[chain] = {}
#                 x = float(line[31:38])
#                 y = float(line[39:46])
#                 z = float(line[47:54])
#                 resnum = int(line[22:26])
#                 atom_name = line[12:16].strip()
#                 # if atom_name in ['CA', 'N', 'C', 'O', 'CB']:
#                 if chain == 'A' and resnum in restraints_a:
#                     if resnum not in coord_dic[chain]:
#                         coord_dic[chain][resnum] = []
#                     coord_dic[chain][resnum].append((x, y, z))
#                 if chain == 'B' and resnum in restraints_b:
#                     if resnum not in coord_dic[chain]:
#                         coord_dic[chain][resnum] = []
#                     coord_dic[chain][resnum].append((x, y, z))
#     coord_a = coord_dic['A']
#     coord_b = coord_dic['B']
#     contacts = {'A': [], 'B': []}
#     for res_a in coord_a:
#         for atom_coords_a in coord_a[res_a]:
#             x_a, y_a, z_a = atom_coords_a
#             for res_b in coord_b:
#                 for atom_coords_b in coord_b[res_b]:
#                     x_b, y_b, z_b = atom_coords_b
#                     c_i = x_a, y_a, z_a
#                     c_j = x_b, y_b, z_b
#                     distance = scipy.spatial.distance.euclidean(c_i, c_j)
#                     if distance <= cutoff and distance > 2.0:
#                         # assumme that if distance is > 2.0 it is a clash
#                         contacts['A'].append(res_a)
#                         contacts['B'].append(res_b)
#     # check how many are satisfied
#     satisfied_a = len(set(contacts['A']).intersection(restraints_a))
#     satisfied_b = len(set(contacts['B']).intersection(restraints_b))
#     total_restraints = len(restraints_a) + len(restraints_b)
#     satisfied_ratio = (satisfied_a + satisfied_b) / total_restraints
#     return satisfied_ratio
