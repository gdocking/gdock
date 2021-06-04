import configparser
import os
import logging
import math
import pathlib
import multiprocessing
import sys
from utils.files import get_full_path
from modules.profit import Profit
etc_folder = get_full_path('etc')
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, 'gdock.ini'), encoding='utf-8')
fcc_path = ini.get('third_party', 'fcc_path')
sys.path.append(f"{fcc_path}/scripts")
ga_log = logging.getLogger('ga_log')

try:
    import calc_fcc_matrix
    import make_contacts
    import cluster_fcc
except Exception as e:
    ga_log.error('FCC could not be imported')
    ga_log.error(e)
    sys.exit()


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
        #: FCC things
        self.contact_executable = f"{fcc_path}/src/contact_fcc"

    @staticmethod
    def get_structures(data_dic):
        structure_score_l = []
        for gen in data_dic:
            for ind in data_dic[gen]:
                struct = data_dic[gen][ind]['structure']
                if struct:
                    if 'score' not in data_dic[gen][ind]:
                        raise Exception('Structure does not contain score,'
                                        ' cannot proceed.')
                    score = data_dic[gen][ind]['score']
                    structure_score_l.append((struct, score))
        sorted_structure_score_l = sorted(structure_score_l,
                                          key=lambda x: x[1])

        # get only the structures
        structure_l = [e[0] for e in sorted_structure_score_l]
        return structure_l

    @staticmethod
    def calc_contact(executable, pdb_f, cutoff=5.0):
        # fcc-wrapper
        contact_f = pdb_f.replace('.pdb', '.contacts')
        make_contacts._calculate_contacts(executable, pdb_f, str(cutoff))
        if os.path.isfile(contact_f):
            return contact_f
        else:
            return ''

    def cluster(self, cutoff=0.75, top=1000):
        """Use FCC to cluster structures."""

        # Calculate contacts
        ga_log.info('FCC - Calculating contacts')
        pool = multiprocessing.Pool(processes=self.nproc)
        input_structure_l = self.structure_list[:top]
        for pdb in input_structure_l:
            pool.apply_async(self.calc_contact, args=(self.contact_executable,
                                                      pdb))

        pool.close()
        pool.join()

        contact_file_l = []
        for pdb in input_structure_l:
            contact_f = pdb.replace('.pdb', '.contacts')
            if os.path.isfile(contact_f):
                contact_file_l.append(contact_f)

        # Calculate matrix
        ga_log.info('FCC - Calculating matrix')
        parsed_contacts = calc_fcc_matrix.parse_contact_file(contact_file_l,
                                                             False)
        matrix = calc_fcc_matrix.calculate_pairwise_matrix(parsed_contacts,
                                                           False)

        # write it to a file, so we can read it afterwards and don't
        #  need to reinvent the wheel
        fcc_matrix_f = f'{self.analysis_path}/fcc.matrix'
        with open(fcc_matrix_f, 'w') as fh:
            for data in list(matrix):
                data_str = f"{data[0]} {data[1]} {data[2]:.2f} {data[3]:.3f}"
                data_str += os.linesep
                fh.write(data_str)
        fh.close()

        # cluster
        ga_log.info('FCC - Clustering')
        pool = cluster_fcc.read_matrix(fcc_matrix_f, cutoff, strictness=0.75)
        element_pool, clusters = cluster_fcc.cluster_elements(pool, 4)

        if clusters:
            ga_log.info(f'FCC - {len(self.cluster_dic)} clusters identified')
            # use fcc's output
            cluster_out = f'{self.analysis_path}/cluster.out'
            with open(cluster_out, 'w') as fh:
                cluster_fcc.output_clusters(fh, clusters)
            fh.close()

            # read it again!
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
        else:
            ga_log.info('FCC - No clusters identified')

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

                    generation_str = str(gen).rjust(4, '0')
                    individual_str = str(ind).rjust(4, '0')

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
