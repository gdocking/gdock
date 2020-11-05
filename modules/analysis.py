import configparser
import os
import numpy as np
from scipy.spatial.transform import Rotation as R
from utils.files import get_full_path
from utils.functions import format_coords

etc_folder = get_full_path('etc')
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, 'gdock.ini'), encoding='utf-8')
fcc = ini.get('third_party', 'dockq_exe')


class Analysis:

    def __init__(self, input_structures, result_dic, run_params):
        """Setup the analysis class."""
        self.input_structures = input_structures
        self.result_dic = result_dic
        self.analysis_path = f"{run_params['folder']}/analysis"

    def generate_structures(self):
        """Read information from simulation and generate structures."""
        structure_list = []
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
                structure_list.append(pdb_name)

        return structure_list

    def cluster(self):
        # Use FCC for this
        pass

    def output(self):
        pass

    def plot(self):
        pass
    #     """
    #
    #     :type plot_name: string
    #     """
    #     import pandas as pd
    #     import seaborn as sns
    #     import matplotlib.pyplot as plt
    #     sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    #
    #     def label(color, plot_label):
    #         ax = plt.gca()
    #         ax.text(0, .2, plot_label, fontweight="bold", color=color,
    #                 ha="left", va="center", transform=ax.transAxes)
    #
    #     result_l = []
    #     for gen in self.generation_dic:
    #         for ind in self.generation_dic[gen]:
    #             name = 'pdbs/gd_' + '_'.join(map("{:.2f}".format, self.generation_dic[gen][ind][0])) + '.pdb'
    #             fitness = self.generation_dic[gen][ind][1][0]
    #             result_l.append((gen, ind, fitness, name))
    #     df = pd.DataFrame(result_l)
    #     df.columns = ['generation', 'individual', 'fitness', 'pdb']
    #
    #     pal = sns.cubehelix_palette(n_colors=len(df.generation.unique()), rot=0, light=.7, reverse=False)
    #     g = sns.FacetGrid(df, row="generation", hue="generation", aspect=15, height=.5, palette=pal)
    #     g.map(sns.kdeplot, "fitness", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
    #     g.map(sns.kdeplot, "fitness", clip_on=False, color="w", lw=2, bw=.2)
    #     g.map(plt.axhline, y=0, lw=2, clip_on=False)
    #
    #     g.map(label, "fitness")
    #     # Set the subplots to overlap
    #     g.fig.subplots_adjust(hspace=-.25)
    #
    #     # Remove axes details that don't play well with overlap
    #     g.set_titles("")
    #     g.set(yticks=[])
    #     g.despine(bottom=True, left=True)
    #     plt.savefig(plot_name)

    # def output(self):
    #     """Output the fitness and individual properties."""
    #     folder = self.run_params['folder']
    #     output_f = f'{folder}/gdock.dat'
    #     ga_log.info(f'Writing output to {output_f}')
    #     with open(output_f, 'w') as fh:
    #         for gen in self.generation_dic:
    #             for ind in self.generation_dic[gen]:
    #                 individual = ' '.join(map(str, self.generation_dic[gen][ind][0]))
    #                 fitness = self.generation_dic[gen][ind][1][0]
    #                 tbw = f'{gen},{ind},{fitness},{individual}\n'
    #                 fh.write(tbw)
    #     fh.close()
