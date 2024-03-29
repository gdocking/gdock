import argparse
import glob
import logging
import os
import shutil
from pathlib import Path

bm_log = logging.getLogger("bm_log")
bm_log.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter(
    " %(asctime)s %(module)s:%(lineno)d" " %(levelname)s - %(message)s"
)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "bm5_path", help=("Location of the Protein-Protein" " Docking Benchmark v5")
    )
    parser.add_argument(
        "gdockbm_path", help=("Location where the prepared" " folders will be")
    )
    parser.add_argument(
        "--np",
        help=("Number of processors to be used for" " each run"),
        type=int,
        default=4,
    )
    args = parser.parse_args()

    folder_list = []
    to_skip = [".", "ana_scripts", "data", "scripts"]
    for folder in glob.glob(f"{args.bm5_path}/HADDOCK-ready/*"):
        folder_name = Path(folder).name
        if not any([j in folder_name for j in to_skip]):
            folder_list.append(folder)

    folder_list.sort()

    bm_log.info(f"Found {len(folder_list)} targets")

    if not os.path.isdir(args.gdockbm_path):
        bm_log.info(f"Creating Benchmark Path {args.gdockbm_path}")
        os.mkdir(args.gdockbm_path)
    else:
        bm_log.info(f"Path already created {args.gdockbm_path}")

    for target in folder_list:
        bm_log.info(f"Preparing {target}")
        loc = Path(target)
        target_path = str(loc.parent)
        target_name = loc.name

        contact_file = f"{target_path}/{target_name}/ana_scripts/" "target.contacts3.9"
        izone_file = f"{target_path}/{target_name}/ana_scripts/target.izone"
        receptor_unbound = f"{target_name}_r_u.pdb"
        ligand_unbound = f"{target_name}_l_u.pdb"
        native = "ana_scripts/target.pdb"

        bm_log.debug(f"contact_file {contact_file}")
        bm_log.debug(f"izone_file {izone_file}")
        bm_log.debug(f"receptor {receptor_unbound}")
        bm_log.debug(f"ligand {ligand_unbound}")
        bm_log.debug(f"native {native}")

        restraints = get_restraints(contact_file)

        # Prepare folder
        if os.path.isdir(f"{args.gdockbm_path}/{target_name}"):
            bm_log.warning(f"{target_name} found in {args.gdockbm_path}")
            shutil.rmtree(f"{args.gdockbm_path}/{target_name}")

        bm_log.info(f"Creating {args.gdockbm_path}/{target_name}")
        os.mkdir(f"{args.gdockbm_path}/{target_name}")

        bm_log.info("Copying needed files to the target folder")
        # receptor
        shutil.copy(
            f"{target_path}/{target_name}/{receptor_unbound}",
            f"{args.gdockbm_path}/{target_name}",
        )
        # ligand
        shutil.copy(
            f"{target_path}/{target_name}/{ligand_unbound}",
            f"{args.gdockbm_path}/{target_name}",
        )
        # native
        shutil.copy(
            f"{target_path}/{target_name}/{native}",
            (f"{args.gdockbm_path}/{target_name}/" f"{target_name}_complex_bound.pdb"),
        )
        # izone
        shutil.copy(
            f"{target_path}/{target_name}/ana_scripts/target.izone",
            (
                f"{args.gdockbm_path}/{target_name}/"
                f"{target_name}_complex_bound.izone"
            ),
        )

        # Write run.toml
        run_toml_f = f"{args.gdockbm_path}/{target_name}/run.toml"
        bm_log.info(f"Writing run parameters to {run_toml_f}")
        run = "[main]" + os.linesep
        run += "identifier = 'run'" + os.linesep
        run += f"number_of_processors = {args.np}" + os.linesep
        run += os.linesep
        run += "[restraints]" + os.linesep
        run += f"A = {restraints['A']}" + os.linesep
        run += f"B = {restraints['B']}" + os.linesep
        run += os.linesep
        run += "[molecules]" + os.linesep
        run += f"A = '{receptor_unbound}'" + os.linesep
        run += f"B = '{ligand_unbound}'" + os.linesep
        run += f"native = '{target_name}_complex_bound.pdb'" + os.linesep

        with open(run_toml_f, "w") as run_fh:
            run_fh.write(run)
        run_fh.close()

        bm_log.info(f"Setup of {target_name} complete")
