# WIP
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# def plot(self):
#     """Create a ridge plot to visually check on the generations"""
#     plot_name = f'{self.analysis_path}/plot.png'
#     sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
#
#     def label(color, plot_label):
#         ax = plt.gca()
#         ax.text(0, .2, plot_label, fontweight="bold", color=color,
#                 ha="left", va="center", transform=ax.transAxes)
#
#     result_l = []
#     for gen in self.result_dic:
#         for ind in self.result_dic[gen]:
#             _, fitness = self.result_dic[gen][ind]
#             fitness = fitness[0]  # placeholder for multiple values of fitnessess
#             name = f"{str(gen).rjust(3, '0')}_{str(ind).rjust(3, '0')}"
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