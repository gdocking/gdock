"""Functions Module."""
import logging
import shlex
import tempfile
import subprocess  # nosec
import os
import sys
import numpy as np

ga_log = logging.getLogger('ga_log')


def get_coords(pdb_f, target_chain=None):
    """Read PDB file and return array with all atoms."""
    coord = []
    with open(pdb_f, 'r') as file_handle:
        for line in file_handle.readlines():
            if line.startswith('ATOM'):
                chain = line[21]

                x_c = float(line[31:38])
                y_c = float(line[39:46])
                z_c = float(line[47:54])

                if target_chain and target_chain == chain:
                    coord.append((x_c, y_c, z_c))
                elif not target_chain:
                    coord.append((x_c, y_c, z_c))
    return np.array(coord)


def add_dummy(pdb_f, output_f, coor_list):
    """Add a dummy atom to a PDB file according to a list of coordinates."""
    new_pdb = []
    with open(pdb_f, 'r') as ref_fh:
        for line in ref_fh.readlines():
            if line.startswith('ATOM'):
                new_pdb.append(line)
    ref_fh.close()

    for dummy_coord in coor_list:
        dum_x = f'{dummy_coord[0]:.3f}'.rjust(7, ' ')
        dum_y = f'{dummy_coord[1]:.3f}'.rjust(7, ' ')
        dum_z = f'{dummy_coord[2]:.3f}'.rjust(7, ' ')
        dummy_line = \
          f'ATOM    999  H   DUM X   1     {dum_x} {dum_y} {dum_z}  1.00  1.00           H  \n'
        new_pdb.append(dummy_line)

    with open(output_f, 'w') as out_fh:
        for line in new_pdb:
            out_fh.write(line)
    out_fh.close()
    return True


def tidy(pdb_str):
    """Save temporary file and retrieve it as string."""
    tmp = tempfile.NamedTemporaryFile()
    tmp_out = tempfile.NamedTemporaryFile()
    with open(tmp.name, 'w') as file_handle:
        file_handle.write(pdb_str)
    cmd = f'pdb_tidy {tmp.name}'
    ga_log.debug('Tidying up with command %s', cmd)
    out = open(f'{tmp_out.name}', 'w')
    proc = subprocess.Popen(shlex.split(cmd),
                            shell=False,
                            stdout=out,
                            stderr=subprocess.PIPE)  # nosec
    proc.communicate()
    if not os.path.isfile(tmp.name):
        ga_log.error('Could not tidy the pdb!')
        sys.exit()
    else:
        tidy_pdb_str = tmp_out.read()
        return tidy_pdb_str.decode()


def write_coords(pdb_f, output_f, coords):
    """Read a PDB and rewrite it using a coordinate list."""
    counter = 0
    with open(output_f, 'w') as out_fh:
        with open(pdb_f, 'r') as ref_fh:
            for line in ref_fh.readlines():
                if line.startswith('ATOM'):
                    new_x = f'{coords[counter][0]:.3f}'.rjust(7, ' ')
                    new_y = f'{coords[counter][1]:.3f}'.rjust(7, ' ')
                    new_z = f'{coords[counter][2]:.3f}'.rjust(7, ' ')
                    new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}'
                    out_fh.write(new_line)
                    counter += 1
        ref_fh.close()
    out_fh.close()
    return True


def draw_dummy(output_f, dummy_coord):
    """Create a dummy atom in space."""
    with open(output_f, 'w') as out_fh:
        dum_x = f'{dummy_coord[0]:.3f}'.rjust(7, ' ')
        dum_y = f'{dummy_coord[1]:.3f}'.rjust(7, ' ')
        dum_z = f'{dummy_coord[2]:.3f}'.rjust(7, ' ')
        dummy_line = \
            f'ATOM    999  H   DUM X   1     {dum_x} {dum_y} {dum_z}     1.00  1.00           H  \n'
        out_fh.write(dummy_line)
    out_fh.close()


def format_coords(coord):
    """Make a set of coordinated PDB-format ready."""
    new_x = f'{coord[0]:.3f}'.rjust(7, ' ')
    new_y = f'{coord[1]:.3f}'.rjust(7, ' ')
    new_z = f'{coord[2]:.3f}'.rjust(7, ' ')
    return new_x, new_y, new_z
