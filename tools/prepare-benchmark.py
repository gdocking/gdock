# Prepare the benchmark
import glob
import os
import shutil
from pathlib import Path
import logging
import argparse

bm_log = logging.getLogger('bm_log')
bm_log.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter(' %(asctime)s %(module)s:%(lineno)d %(levelname)s - %(message)s')
ch.setFormatter(formatter)
bm_log.addHandler(ch)


def get_restraints(contact_f):
    """Parse a contact file and retrieve the restraints."""
    restraint_dic = {}
    with open(contact_f) as fh:
        for line in fh.readlines():
            if not line:  # could be empty
                continue
            data = line.split()
            res_i = int(data[0])
            chain_i = data[1]
            res_j = int(data[3])
            chain_j = data[4]

            if chain_i not in restraint_dic:
                restraint_dic[chain_i] = []
            if chain_j not in restraint_dic:
                restraint_dic[chain_j] = []

            if res_i not in restraint_dic[chain_i]:
                restraint_dic[chain_i].append(res_i)
            if res_j not in restraint_dic[chain_j]:
                restraint_dic[chain_j].append(res_j)

    return restraint_dic


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("bm5_path", help="Location of the Protein-Protein Docking Benchmark v5")
    parser.add_argument("gdockbm_path", help="Location where the prepared folders will be")
    args = parser.parse_args()

    target_folders = glob.glob(f'{args.bm5_path}/HADDOCK-ready/*')
    excluded_folders = ['ana_scripts', 'data', 'scripts']

    folder_list = list(set(target_folders) - set(excluded_folders))
    bm_log.info(f'Found {len(folder_list)} targets')

    if not os.path.isdir(args.gdockbm_path):
        bm_log.info(f'Creating Benchmark Path {args.gdockbm_path}')
        os.path.mkdir(args.gdockbm_path)
    else:
        bm_log.info(f'Path already created {args.gdockbm_path}')

    for target in folder_list:
        bm_log.info(f'Preparing {target}')
        loc = Path(target)
        target_path = str(loc.parent)
        target_name = loc.name

        contact_file = f'{target_path}/{target_name}/ana_scripts/target.contacts3.9'
        receptor_unbound = f'{target_name}_r_u.pdb'
        ligand_unbound = f'{target_name}_l_u.pdb'
        native = 'ana_scripts/target.pdb'

        bm_log.debug(f'contact_file {contact_file}')
        bm_log.debug(f'receptor {receptor_unbound}')
        bm_log.debug(f'ligand {ligand_unbound}')
        bm_log.debug(f'native {native}')

        restraints = get_restraints(contact_file)

        # Prepare folder
        if os.path.isdir(f'{args.gdockbm_path}/{target_name}'):
            bm_log.warning(f'{target_name} found in {args.gdockbm_path}')
            # shutil.rmtree(f'{args.gdockbm_path}/{target_name}')
        else:
            bm_log.info(f'Creating {args.gdockbm_path}/{target_name}')
            os.mkdir(f'{args.gdockbm_path}/{target_name}')

        bm_log.info('Copying receptor, ligand and native files to the target folder')
        shutil.copy(f'{target_path}/{target_name}/{receptor_unbound}', f'{args.gdockbm_path}/{target_name}')
        shutil.copy(f'{target_path}/{target_name}/{ligand_unbound}', f'{args.gdockbm_path}/{target_name}')
        shutil.copy(f'{target_path}/{target_name}/{native}', f'{args.gdockbm_path}/{target_name}/native.pdb')

        # Write run.toml
        bm_log.info(f'Writing run parameters to {args.gdockbm_path}/{target_name}/run.toml')
        run = "[main]\n"
        run += "identifier = 'run'\n"
        run += "number_of_processors = 4\n"
        run += "\n"
        run += "[restraints]\n"
        run += f"A = {restraints['A']}\n"
        run += f"B = {restraints['B']}\n"
        run += "\n"
        run += "[molecules]\n"
        run += f"A = '{receptor_unbound}'\n"
        run += f"B = '{ligand_unbound}'\n"
        run += "native = 'native.pdb'"

        with open(f'{args.gdockbm_path}/{target_name}/run.toml', 'w') as run_fh:
            run_fh.write(run)
        run_fh.close()

        bm_log.info(f'Setup of {target_name} complete')
