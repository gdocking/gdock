# from a complex, extract its true-interface
import argparse
import configparser
import os
import subprocess
from utils.files import get_full_path

etc_folder = get_full_path('etc')
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, 'gadock.ini'), encoding='utf-8')
contact_exe = ini.get('third_party', 'contact_exe')


def get_ti(pdb_f, cutoff=4.9):
    """Get the True-Interface of a PDB complex

    :param pdb_f:
    :param cutoff:
    :return:
    """
    ti_dic = {}
    cmd = f'{contact_exe} {pdb_f} {cutoff}'
    proc = subprocess.Popen(cmd.split(), shell=False, stdout=subprocess.PIPE)
    out = proc.stdout.read().decode('utf-8').split('\n')
    for line in out:
        if not line:  # could be empty
            continue
        data = line.split()
        res_i = int(data[0])
        chain_i = data[1]
        res_j = int(data[3])
        chain_j = data[4]

        if chain_i not in ti_dic:
            ti_dic[chain_i] = []
        if chain_j not in ti_dic:
            ti_dic[chain_j] = []

        if res_i not in ti_dic[chain_i]:
            ti_dic[chain_i].append(res_i)
        if res_j not in ti_dic[chain_j]:
            ti_dic[chain_j].append(res_j)

    return ti_dic


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_pdb')

    args = parser.parse_args()

    interface_dic = get_ti(args.input_pdb)

    for c in interface_dic:
        reslist = interface_dic[c]
        reslist.sort()
        reslist_str = ','.join(map(str, reslist))
        print(f'{c} = [{reslist_str}]')\
