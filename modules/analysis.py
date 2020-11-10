import configparser
import os
import subprocess  # nosec
import shlex
import logging
import numpy as np
from scipy.spatial.transform import Rotation as R
from utils.files import get_full_path
from utils.functions import format_coords

ga_log = logging.getLogger('ga_log')
etc_folder = get_full_path('etc')
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, 'gdock.ini'), encoding='utf-8')
fcc = ini.get('third_party', 'fcc_path')


class Analysis:

    def __init__(self, input_structures, result_dic, run_params):
        """Setup the analysis class."""
        self.input_structures = input_structures
        self.result_dic = result_dic
        self.analysis_path = f"{run_params['folder']}/analysis"
        self.nproc = run_params['np']
        self.structure_list = []
        self.cluster_dic = {}

    def generate_structures(self):
        """Read information from simulation and generate structures."""
        input_structure_dic = {}
        for line in self.input_structures.split('\n'):
            if line.startswith('ATOM'):
                chain = line[21]
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                if chain not in input_structure_dic:
                    input_structure_dic[chain] = {'coord': [], 'raw': []}
                input_structure_dic[chain]['coord'].append((x, y, z))
                input_structure_dic[chain]['raw'].append(line)

        # use the chromossome and create the structure!
        for gen in self.result_dic:
            for idx in self.result_dic[gen]:
                individual = self.result_dic[gen][idx][0]

                # transform
                rot = R.from_euler('zyx', individual[:3])
                transl = individual[3:]
                c = np.array(input_structure_dic['B']['coord'])
                center = c.mean(axis=0)
                c -= center
                c -= transl
                r = np.array([rot.apply(e) for e in c])
                r += center
                input_structure_dic['B']['coord'] = list(r)

                pdb_name = f"{self.analysis_path}/{str(gen).rjust(3, '0')}_{str(idx).rjust(3, '0')}.pdb"
                with open(pdb_name, 'w') as fh:
                    for chain in input_structure_dic:
                        for coord, line in zip(input_structure_dic[chain]['coord'], input_structure_dic[chain]['raw']):
                            new_x, new_y, new_z = format_coords(coord)
                            new_line = f'{line[:30]} {new_x} {new_y} {new_z} {line[55:]}\n'
                            fh.write(new_line)
                fh.close()
                self.structure_list.append(pdb_name)

        self.structure_list.sort()
        return self.structure_list

    def cluster(self, cutoff=0.75):
        """Use FCC to cluster structures."""
        ga_log.info('FCC - making contacts')
        # TODO: Make this run using multiple processors
        make_contacts = f'{fcc}/scripts/make_contacts.py'
        pdb_list = f'{self.analysis_path}/pdb.list'
        with open(pdb_list, 'w') as fh:
            for pdb in self.structure_list:
                fh.write(f'{pdb}\n')
        fh.close()

        cmd = f'{make_contacts} -f {pdb_list}'
        ga_log.debug(f'cmd is: {cmd}')
        out = subprocess.check_output(shlex.split(cmd), shell=False, stderr=subprocess.PIPE)  # nosec
        if 'Finished' not in out.decode('utf-8'):
            ga_log.error('FCC - make_contacts.py failed')
            exit()

        ga_log.info('FCC - Calculating matrix')
        calc_fcc_matrix = f'{fcc}/scripts/calc_fcc_matrix.py'
        contact_list = f'{self.analysis_path}/contact.list'
        fcc_matrix = f'{self.analysis_path}/fcc.matrix'
        with open(contact_list, 'w') as fh:
            for pdb in self.structure_list:
                contact_fname = pdb.replace('.pdb', '.contacts')
                fh.write(f'{contact_fname}\n')
        fh.close()

        cmd = f'{calc_fcc_matrix} -f {contact_list} -o {fcc_matrix}'
        subprocess.check_output(shlex.split(cmd), shell=False, stderr=subprocess.PIPE)  # nosec
        if not os.path.isfile(fcc_matrix):
            ga_log.error('FCC - calc_fcc_matrix.py failed')
            exit()

        ga_log.info(f'FCC - Clustering with cutoff={cutoff}')
        cluster_fcc = f'{fcc}/scripts/cluster_fcc.py'
        cluster_out = f'{self.analysis_path}/cluster.out'
        cmd = f'{cluster_fcc} {fcc_matrix} {cutoff} -o {cluster_out}'
        subprocess.check_output(shlex.split(cmd), shell=False, stderr=subprocess.PIPE)  # nosec
        if os.stat(cluster_out).st_size == 0:
            ga_log.warning('No clusters were written!')
            return
        else:
            with open(cluster_out, 'r') as fh:
                for line in fh.readlines():
                    data = line.split()
                    cluster_id = int(data[1])
                    self.cluster_dic[cluster_id] = []
                    cluster_elements = map(int, data[4:])
                    for element in cluster_elements:
                        element_name = os.path.split(self.structure_list[element])[1].split('.pdb')[0]
                        self.cluster_dic[cluster_id].append(element_name)
            ga_log.info(f'FCC - {len(self.cluster_dic)} clusters identified')

        return self.cluster_dic

    def output(self):
        """Generate a script friendly output table."""
        output_f = f'{self.analysis_path}/gdock.dat'
        ga_log.info(f'Saving output file to {output_f}')
        with open(output_f, 'w') as fh:
            fh.write('generation,individual,fitness\n')
            for generation in self.result_dic:
                for individual in self.result_dic[generation]:
                    _, fitness = self.result_dic[generation][individual]
                    fitness = fitness[0]  # placeholder for multiple values of fitnessess
                    fh.write(f"{str(generation).rjust(3, '0')},{str(individual).rjust(3, '0')},{fitness}\n")
        fh.close()
