import configparser
import os
import sys
import subprocess  # nosec
import shlex
import logging
import math
import pathlib
from utils.files import get_full_path
from modules.profit import Profit

ga_log = logging.getLogger('ga_log')
etc_folder = get_full_path('etc')
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, 'gdock.ini'), encoding='utf-8')
fcc = ini.get('third_party', 'fcc_path')


class Analysis:

    def __init__(self, result_dic, run_params):
        """Setup the analysis class."""
        self.result_dic = result_dic
        self.analysis_path = f"{run_params['folder']}/analysis"
        if 'native' in run_params:
            self.native = run_params['native']
        else:
            self.native = ''
        self.nproc = run_params['np']
        self.structure_list = self.get_structures(result_dic)
        self.cluster_dic = {}
        self.irmsd_dic = {}

    @staticmethod
    def get_structures(data_dic):
        structure_l = []
        for gen in data_dic:
            for ind in data_dic[gen]:
                struct = data_dic[gen][ind]['structure']
                if struct:
                    structure_l.append(struct)
        return structure_l

    def cluster(self, cutoff=0.75):
        """Use FCC to cluster structures."""
        ga_log.info('FCC - Calculating contacts')
        make_contacts = f'{fcc}/scripts/make_contacts.py'
        pdb_list = f'{self.analysis_path}/pdb.list'
        with open(pdb_list, 'w') as fh:
            for pdb in self.structure_list:
                fh.write(f'{pdb}' + os.linesep)
        fh.close()

        cmd = f'{make_contacts} -n {self.nproc}  -f {pdb_list}'
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
            ga_log.warning('No clusters were found')
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

        sep = '\t'
        header = 'gen' + sep
        header += 'ind' + sep
        header += 'ranking' + sep
        header += 'score' + sep
        header += 'fitness' + sep
        header += 'energy' + sep
        header += 'irmsd' + sep
        header += 'cluster_id' + sep
        header += 'internal_cluster_ranking' + sep
        header += 'structure_path' + os.linesep

        with open(output_f, 'w') as fh:
            fh.write(header)
            for gen in self.result_dic:
                for ind in self.result_dic[gen]:
                    fitness_values = self.result_dic[gen][ind]['fitness']

                    # placeholder for multiple values of fitnessess
                    fitness = fitness_values[0]

                    score = self.result_dic[gen][ind]['score']
                    ranking = self.result_dic[gen][ind]['ranking']
                    energy = self.result_dic[gen][ind]['energy']

                    generation_str = str(gen).rjust(3, '0')
                    individual_str = str(ind).rjust(3, '0')

                    structure = self.result_dic[gen][ind]['structure']
                    if structure:
                        structure_id = pathlib.Path(structure).stem
                        try:
                            irmsd = self.irmsd_dic[structure_id]
                        except KeyError:
                            irmsd = float('nan')
                    else:
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

                    output_str = f"{gen}" + sep
                    output_str += f"{ind}" + sep
                    output_str += f"{ranking}" + sep
                    output_str += f"{score:.3f}" + sep
                    output_str += f"{fitness:.3f}" + sep
                    output_str += f"{energy:.3f}" + sep
                    output_str += f"{irmsd:.2f}" + sep
                    output_str += f"{cluster_id}" + sep
                    output_str += f"{internal_ranking}" + sep
                    output_str += f"{structure}" + os.linesep
                    ga_log.debug(output_str)

                    if not math.isnan(ranking):
                        # do not write to data file if there's no ranking
                        fh.write(output_str)

        fh.close()
