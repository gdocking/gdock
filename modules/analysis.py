
class Analysis:

    def __init__(self, input_structures, result_dic):
        """Setup the analysis class."""
        self.input_structures = input_structures
        self.result_dic = result_dic

    def generate_structures(self):
        """Read information from simulation and generate structures"""
        pass

    def cluster(self):
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
