import ast
import configparser
import logging
import os
import secrets
import shlex
import subprocess  # nosec
import sys
import tempfile
from pathlib import Path

import numpy

from utils.files import get_full_path

ga_log = logging.getLogger("ga_log")

etc_folder = get_full_path("etc")
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, "gdock.ini"), encoding="utf-8")
pdbtools_path = ini.get("third_party", "pdbtools_path")
python_exe = sys.executable


def tidy(pdb_str):
    """Tidy PDB using pdbtools pdb_tidy."""
    input_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    input_pdb.write(str.encode(pdb_str))
    input_pdb.close()

    cmd = f"{python_exe} {pdbtools_path}/pdbtools/pdb_tidy.py {input_pdb.name}"
    ga_log.debug(f"Tidying up with command {cmd}")

    output = subprocess.check_output(shlex.split(cmd))  # nosec

    os.unlink(input_pdb.name)

    if not output:
        ga_log.error("Could not tidy the pdb!")
        exit()
    else:
        tidy_pdb = output.decode("utf-8")
        return tidy_pdb


def format_coords(coord):
    """Make a set of coordinated PDB-format ready."""
    new_x = f"{coord[0]:.3f}".rjust(7, " ")
    new_y = f"{coord[1]:.3f}".rjust(7, " ")
    new_z = f"{coord[2]:.3f}".rjust(7, " ")
    return new_x, new_y, new_z


def random_quote():
    """Retrieve a correctly formatted quote."""
    quote_file = f"{get_full_path('etc')}/quotes"
    quote_list = []
    if not os.path.isfile(quote_file):
        return "", ""
    else:
        with open(quote_file, "r") as quote_h:
            for line in quote_h.readlines():
                if not line.startswith("#"):
                    try:
                        auth, quote = line.split("__")
                        quote = quote[:-1]
                        quote_list.append((auth, quote))
                    except ValueError:
                        auth = ""
                        quote = ""

        pick = secrets.choice(range(0, len(quote_list)))
        random_author, random_quote = quote_list[pick]
        return random_author, random_quote


def du(path):
    """Disk Usage."""
    size_in_bytes = sum(file.stat().st_size for file in Path(path).rglob("*"))
    size_in_kb = size_in_bytes / 1024
    size_in_mb = size_in_bytes / (1024 * 1024)
    size_in_gb = size_in_bytes / (1024 * 1024 * 1024)
    if int(size_in_gb) > 0:
        return f"{size_in_gb:.2f} GB"
    elif int(size_in_mb) > 0:
        return f"{size_in_mb:.2f} MB"
    elif int(size_in_kb) > 0:
        return f"{size_in_kb:.2f} KB"
    else:
        return f"{size_in_bytes} B"


def check_if_py3(code_f):
    """Test if the code is python3 compatible."""
    try:
        with open(code_f, "rb") as code_fh:
            ast.parse(code_fh.read())
        code_fh.close()
    except SyntaxError as e:
        ga_log.debug(e)
        return False

    return True


def summary(value_list):
    """Calculate the summary statistics of a value list."""
    mean = numpy.mean(value_list)
    std = numpy.std(value_list)
    max_v = max(value_list)
    min_v = min(value_list)
    return {"mean": mean, "std": std, "max": max_v, "min": min_v}


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
#        dummy_line = (f'ATOM    999  H   DUM X   1     {dum_x} {dum_y}'
#                      f' {dum_z}  1.00  1.00           H  \n')
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
#                    new_line = (f'{line[:30]} {new_x} {new_y} {new_z}'
#                                f' {line[55:]}')
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
#        dummy_line = (f'ATOM    999  H   DUM X   1     {dum_x} {dum_y} '
#                      f'{dum_z}     1.00  1.00           H  \n')
#         out_fh.write(dummy_line)
#     out_fh.close()
