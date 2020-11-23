import configparser
import os
import subprocess  # nosec
import shlex
import logging
import multiprocessing
from pathlib import Path
import numpy as np
from utils.files import get_full_path
from utils.functions import format_coords
from modules.geometry import Geometry

ga_log = logging.getLogger('ga_log')
etc_folder = get_full_path('etc')
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, 'gdock.ini'), encoding='utf-8')
fcc = ini.get('third_party', 'fcc_path')
dockq_exe = ini.get('third_party', 'dockq_exe')


class Analysis:

    def __init__(self, input_structures, result_dic, run_params):
        """Setup the analysis class."""
        self.input_structures = input_structures
        self.result_dic = result_dic
        self.analysis_path = f"{run_params['folder']}/analysis"
        if 'native' in run_params:
            self.native = run_params['native']
        else:
            self.native = ''
        self.nproc = run_params['np']
        self.structure_list = []
        self.cluster_dic = {}
        self.irmsd_dic = {}

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
        c = np.array(input_structure_dic['B']['coord'])
        for gen in self.result_dic:
            for idx in self.result_dic[gen]:

                individual = self.result_dic[gen][idx][0]

                translation_center = individual[3:]
                rotation_angles = individual[:3]

                translated_coords = Geometry.translate(c, translation_center)
                rotated_coords = Geometry.rotate(translated_coords, rotation_angles)

                input_structure_dic['B']['coord'] = list(rotated_coords)

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
        ga_log.info('Calculating contacts')
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

        ga_log.info('Calculating contact matrix')
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
                    cluster_elements = list(map(int, data[4:]))
                    for element in cluster_elements:
                        element_name = os.path.split(self.structure_list[element])[1].split('.pdb')[0]
                        self.cluster_dic[cluster_id].append(element_name)
            ga_log.info(f'FCC - {len(self.cluster_dic)} clusters identified')

        return self.cluster_dic

    def capri_eval(self):
        """Calculate the CAPRI metrics for benchmarking purposes."""
        if not self.native:
            ga_log.warning('Native not defined, capri evaluation skipped')
            return
        ga_log.info(f'Calculating iRMSD against native structure {self.native}')
        pool = multiprocessing.Pool(processes=self.nproc)  # no logging inside the pool because its way to complicated
        results = []
        for pdb in self.structure_list:
            results.append((pdb, pool.apply_async(self.calc_irmsd, args=(dockq_exe, pdb, self.native))))
        for p in results:
            pdb, irmsd = p[0], p[1].get()
            pdb_identifier = Path(pdb).name[:-4]
            self.irmsd_dic[pdb_identifier] = irmsd
        pool.close()
        pool.join()

    def calc_irmsd(self, dockq_exe, pdb_f, native):
        """Calculate the interface root mean square deviation."""
        cmd = f'{dockq_exe} {pdb_f} {native}'
        ga_log.debug(f'cmd is: {cmd}')
        try:
            out = subprocess.check_output(shlex.split(cmd), shell=False)  # nosec
            result = out.decode('utf-8')
            irmsd = float(result.split('\n')[-6].split()[-1])
        except Exception:
            # This will happen when there is something very wrong with the PDB
            irmsd = float('nan')
        return irmsd

    def output(self):
        """Generate a script friendly output table."""
        output_f = f'{self.analysis_path}/gdock.dat'
        ga_log.info(f'Saving output file to {output_f}')
        with open(output_f, 'w') as fh:
            fh.write('generation,individual,fitness,irmsd\n')
            for generation in self.result_dic:
                for individual in self.result_dic[generation]:
                    _, fitness = self.result_dic[generation][individual]
                    fitness = fitness[0]  # placeholder for multiple values of fitnessess
                    generation_str = str(generation).rjust(3, '0')
                    individual_str = str(individual).rjust(3, '0')
                    try:
                        irmsd = self.irmsd_dic[f'{generation_str}_{individual_str}']
                    except KeyError:
                        irmsd = float('nan')
                    fh.write(f"{generation_str},{individual_str},{fitness:.3f},{irmsd:.2f}\n")
        fh.close()
