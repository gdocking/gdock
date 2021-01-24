import shlex
import tempfile
import subprocess  # nosec
import os
import logging
import secrets
from pathlib import Path
from utils.files import get_full_path

ga_log = logging.getLogger('ga_log')


def tidy(pdb_str):
    """Save temporary file and retrieve it as string."""
    tmp = tempfile.NamedTemporaryFile()
    tmp_out = tempfile.NamedTemporaryFile()
    with open(tmp.name, 'w') as f:
        f.write(pdb_str)
    cmd = f'pdb_tidy {tmp.name}'
    ga_log.debug(f'Tidying up with command {cmd}')
    out = open(f'{tmp_out.name}', 'w')
    p = subprocess.Popen(shlex.split(cmd), shell=False, stdout=out, stderr=subprocess.PIPE)  # nosec
    p.communicate()
    if not os.path.isfile(tmp.name):
        ga_log.error('Could not tidy the pdb!')
        exit()
    else:
        tidy_pdb_str = tmp_out.read()
        return tidy_pdb_str.decode()


def format_coords(coord):
    """Make a set of coordinated PDB-format ready."""
    new_x = f'{coord[0]:.3f}'.rjust(7, ' ')
    new_y = f'{coord[1]:.3f}'.rjust(7, ' ')
    new_z = f'{coord[2]:.3f}'.rjust(7, ' ')
    return new_x, new_y, new_z


def random_quote():
    """Retrieve a correctly formatted quote."""
    quote_file = f"{get_full_path('etc')}/quotes"
    quote_list = []
    if not os.path.isfile(quote_file):
        return '', ''
    else:
        with open(quote_file, 'r') as quote_h:
            for line in quote_h.readlines():
                auth, quote = line.split('__')
                quote = quote[:-1]
                quote_list.append((auth, quote))

        random_author, random_quote = quote_list[secrets.choice(range(0, len(quote_list)))]
        return random_author, random_quote


def du(path):
    """Disk Usage."""
    return sum(file.stat().st_size for file in Path(path).rglob('*'))


# ======
#  Helper functions to assist in geometric stuff, dev only
#   Keep it here for prosperity
#
# import numpy as np
#
# def add_dummy(pdb_f, output_f, coor_list):
#     """Add a dummy atom to a PDB file according to a list of coordinates."""
#     new_pdb = []
#     with open(pdb_f, 'r') as ref_fh:
#         for line in ref_fh.readlines():
#             if line.startswith('ATOM'):
#                 new_pdb.append(line)
#     ref_fh.close()
#     for dummy_coord in coor_list:
#         dum_x = f'{dummy_coord[0]:.3f}'.rjust(7, ' ')
#         dum_y = f'{dummy_coord[1]:.3f}'.rjust(7, ' ')
#         dum_z = f'{dummy_coord[2]:.3f}'.rjust(7, ' ')
#         dummy_line = f'ATOM    999  H   DUM X   1     {dum_x} {dum_y} {dum_z}  1.00  1.00           H  \n'
#         new_pdb.append(dummy_line)
#     with open(output_f, 'w') as out_fh:
#         for line in new_pdb:
#             out_fh.write(line)
#     out_fh.close()
#     return True
#
# def write_coords(pdb_f, output_f, coords):
#     """Read a PDB and rewrite it using a coordinate list."""
#     c = 0
#     with open(output_f, 'w') as out_fh:
#         with open(pdb_f, 'r') as ref_fh:
#             for line in ref_fh.readlines():
#                 if line.startswith('ATOM'):
#                     new_x = f'{coords[c][0]:.3f}'.rjust(7, ' ')
#                     new_y = f'{coords[c][1]:.3f}'.rjust(7, ' ')
#                     new_z = f'{coords[c][2]:.3f}'.rjust(7, ' ')
#                     new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}'
#                     out_fh.write(new_line)
#                     c += 1
#         ref_fh.close()
#     out_fh.close()
#     return True
#
# def draw_dummy(output_f, dummy_coord):
#     """Create a dummy atom in space."""
#     with open(output_f, 'w') as out_fh:
#         dum_x = f'{dummy_coord[0]:.3f}'.rjust(7, ' ')
#         dum_y = f'{dummy_coord[1]:.3f}'.rjust(7, ' ')
#         dum_z = f'{dummy_coord[2]:.3f}'.rjust(7, ' ')
#         dummy_line = f'ATOM    999  H   DUM X   1     {dum_x} {dum_y} {dum_z}     1.00  1.00           H  \n'
#         out_fh.write(dummy_line)
#     out_fh.close()
