import configparser
import os
import sys
import subprocess  # nosec
import shlex
import logging
import math
import pathlib
import multiprocessing
import numpy
from utils.files import get_full_path
from utils.functions import format_coords
from modules.geometry import Geometry
from modules.profit import Profit

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
        ga_log.info('Re-generating structures')
        input_structure_dic = {}
        for line in self.input_structures.split(os.linesep):
            if line.startswith('ATOM'):
                chain = line[21]
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                if chain not in input_structure_dic:
                    input_structure_dic[chain] = {'coord': [], 'raw': []}
                input_structure_dic[chain]['coord'].append((x, y, z))
                input_structure_dic[chain]['raw'].append(line)

        pool = multiprocessing.Pool(processes=self.nproc)

        for gen in self.result_dic:
            for idx in self.result_dic[gen]:
                individual = self.result_dic[gen][idx][0]
                pdb_name = (f"{self.analysis_path}/{str(gen).rjust(3, '0')}"
                            f"_{str(idx).rjust(3, '0')}.pdb")
                self.structure_list.append(pdb_name)

                pool.apply_async(self._recreate, args=(input_structure_dic,
                                                       individual,
                                                       pdb_name))

        pool.close()
        pool.join()

    @staticmethod
    def _recreate(input_structure_dic, individual, pdb_name):
        """Use the chromossome information and create the structure."""
        c = numpy.array(input_structure_dic['B']['coord'])

        translation_center = individual[3:]
        rotation_angles = individual[:3]

        translated_coords = Geometry.translate(c, translation_center)
        rotated_coords = Geometry.rotate(translated_coords, rotation_angles)

        input_structure_dic['B']['coord'] = list(rotated_coords)

        with open(pdb_name, 'w') as fh:
            for chain in input_structure_dic:
                coord_l = input_structure_dic[chain]['coord']
                raw_l = input_structure_dic[chain]['raw']
                for coord, line in zip(coord_l, raw_l):
                    new_x, new_y, new_z = format_coords(coord)
                    new_line = (f'{line[:30]} {new_x} {new_y} {new_z}'
                                f' {line[55:]}' + os.linesep)
                    fh.write(new_line)
        fh.close()

    def cluster(self, cutoff=0.75):
        """Use FCC to cluster structures."""
        ga_log.info('FCC - Calculating contacts')
        # TODO: Make this run using multiple processors
        make_contacts = f'{fcc}/scripts/make_contacts.py'
        pdb_list = f'{self.analysis_path}/pdb.list'
        with open(pdb_list, 'w') as fh:
            for pdb in self.structure_list:
                fh.write(f'{pdb}' + os.linesep)
        fh.close()

        cmd = f'{make_contacts} -f {pdb_list}'
        ga_log.debug(f'cmd is: {cmd}')
        try:
            out = subprocess.check_output(shlex.split(cmd),
                                          shell=False,
                                          stderr=subprocess.PIPE)  # nosec
        except subprocess.CalledProcessError as e:
            ga_log.error(f'FCC failed with {e}')
            return
        if 'Finished' not in out.decode('utf-8'):
            ga_log.error('FCC - make_contacts.py failed')
            return

        ga_log.info('FCC - Calculating contact matrix')
        calc_fcc_matrix = f'{fcc}/scripts/calc_fcc_matrix.py'
        contact_list = f'{self.analysis_path}/contact.list'
        fcc_matrix = f'{self.analysis_path}/fcc.matrix'
        with open(contact_list, 'w') as fh:
            for pdb in self.structure_list:
                contact_fname = pdb.replace('.pdb', '.contacts')
                fh.write(f'{contact_fname}' + os.linesep)
        fh.close()

        cmd = f'{calc_fcc_matrix} -f {contact_list} -o {fcc_matrix}'
        subprocess.check_output(shlex.split(cmd),
                                shell=False,
                                stderr=subprocess.PIPE)  # nosec
        if not os.path.isfile(fcc_matrix):
            ga_log.error('FCC - calc_fcc_matrix.py failed')
            sys.exit()

        ga_log.info(f'FCC - Clustering with cutoff={cutoff}')
        cluster_fcc = f'{fcc}/scripts/cluster_fcc.py'
        cluster_out = f'{self.analysis_path}/cluster.out'
        cmd = f'{cluster_fcc} {fcc_matrix} {cutoff} -o {cluster_out}'
        subprocess.check_output(shlex.split(cmd),
                                shell=False,
                                stderr=subprocess.PIPE)  # nosec
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
                        structure_name = self.structure_list[element - 1]
                        element_name = pathlib.Path(structure_name).stem
                        self.cluster_dic[cluster_id].append(element_name)
            ga_log.info(f'FCC - {len(self.cluster_dic)} clusters identified')

        return self.cluster_dic

    def evaluate(self):
        """Use PROFIT to calculate the interface rmsd."""
        try:
            ga_log.info('Using PROFIT to evaluate the irmsd')
            profit = Profit(ref=self.native,
                            mobi=self.structure_list,
                            nproc=self.nproc)
            self.irmsd_dic = profit.calc_irmsd()
        except Exception as e:
            ga_log.warning(e)
            ga_log.warning('Skipping evaluation')

    def output(self):
        """Generate a script friendly output table."""
        output_f = f'{self.analysis_path}/gdock.dat'
        ga_log.info(f'Saving output file to {output_f}')
        with open(output_f, 'w') as fh:
            fh.write('generation,individual,satisfaction,energy,'
                     'irmsd,cluster_id,internal_ranking' + os.linesep)
            for generation in self.result_dic:
                for individual in self.result_dic[generation]:
                    _, fitness = self.result_dic[generation][individual]
                    # placeholder for multiple values of fitnessess
                    satisfaction = fitness[0]
                    energy = fitness[1]
                    generation_str = str(generation).rjust(3, '0')
                    individual_str = str(individual).rjust(3, '0')
                    structure_name = f'{generation_str}_{individual_str}'
                    try:
                        irmsd = self.irmsd_dic[structure_name]
                    except KeyError:
                        irmsd = float('nan')

                    # Add the clustering information
                    # TODO: Improve this part, its very messy.
                    model = f'{generation_str}_{individual_str}'
                    internal_ranking = float('nan')
                    cluster_id = float('nan')
                    if self.cluster_dic:
                        for cluster_id in self.cluster_dic:
                            subcluster = self.cluster_dic[cluster_id]
                            if model in subcluster:
                                # the center has already been removed
                                internal_ranking = subcluster.index(model) + 1
                                break
                    if math.isnan(internal_ranking):
                        cluster_id = float('nan')

                    output_str = f"{generation},"
                    output_str += f"{individual},"
                    output_str += f"{satisfaction:.3f},"
                    output_str += f"{energy:.3f},"
                    output_str += f"{irmsd:.2f},"
                    output_str += f"{cluster_id},"
                    output_str += f"{internal_ranking}" + os.linesep
                    fh.write(output_str)

                    ga_log.debug(output_str)
        fh.close()
