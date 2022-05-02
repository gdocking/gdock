import collections
import logging
import math
import os
import pathlib
import random
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from fccpy import read_pdb
from fccpy.clustering import disjoint_taylor_butina as DTB_ALGO
from fccpy.contacts import (
    BY_RESIDUE,
    get_intermolecular_contacts,
    hash_many,
    read_contacts,
    write_contacts,
)
from fccpy.similarity import build_matrix
from fccpy.similarity import fcc as FCC_METRIC
from fccpy.similarity import write_matrix
from gdock.modules.profit import Profit

matplotlib.use("Agg")
plt.rcParams["font.family"] = "serif"


ga_log = logging.getLogger("ga_log")


class Analysis:
    def __init__(self, result_dic, params):
        """Setup the analysis class."""
        self.result_dic = result_dic
        self.analysis_path = f"{params['folder']}/analysis"
        if "native" in params:
            self.native = params["native"]
        else:
            self.native = ""
        self.nproc = params["main"]["number_of_processors"]
        self.structure_list = self.get_structures(result_dic)
        self.clust_cutoff = params["analysis"]["clust_cutoff"]
        self.clust_min_size = params["analysis"]["clust_min_size"]
        self.clust_n = params["analysis"]["clust_n"]
        self.contact_cutoff = params["analysis"]["contact_cutoff"]
        self.cluster_dic = {}
        self.irmsd_dic = {}

        self._rank_by_energy()

    def _rank_by_energy(self):
        """Add a ranking key to the result dictionary based on the energy."""
        ranked_list = []
        for gen in self.result_dic:
            for ind in self.result_dic[gen]:
                pdb_f = self.result_dic[gen][ind]["structure"]
                if pdb_f:
                    energy = self.result_dic[gen][ind]["energy"]
                    pdb_name = pathlib.Path(pdb_f).stem
                    ranked_list.append((pdb_name, energy))

        ranked = {}
        for counter, (pdb_id, score) in enumerate(
            sorted(ranked_list, key=lambda tup: tup[1]), start=1
        ):
            ranked[pdb_id] = {"rank": counter, "score": score}

        # add the ranking to the data structure
        for gen in self.result_dic:
            for ind in self.result_dic[gen]:

                pdb_f = self.result_dic[gen][ind]["structure"]

                if pdb_f:
                    pdb_id = pathlib.Path(pdb_f).stem
                    ranking = ranked[pdb_id]["rank"]
                else:
                    # there is no structure, its a clone
                    ranking = float("nan")

                self.result_dic[gen][ind]["ranking"] = ranking

    @staticmethod
    def get_structures(data_dic):
        """Retrieve the structures."""
        structure_l = []
        for gen in data_dic:
            for ind in data_dic[gen]:
                struct = data_dic[gen][ind]["structure"]
                if struct:
                    if "energy" not in data_dic[gen][ind]:
                        raise Exception(
                            "Structure does not contain energy, cannot proceed."
                        )
                    satisfaction = data_dic[gen][ind]["satisfaction"]
                    energy = data_dic[gen][ind]["energy"]

                    structure_l.append((struct, satisfaction, energy))

        sorted_structure_l = sorted(
            structure_l, key=lambda x: (x[1], x[2]), reverse=True
        )

        # get only the structures
        structure_l = [e[0] for e in sorted_structure_l]
        return structure_l

    def calc_contact(self, pdb_f):
        """Calculate the intra-molecular contacts."""
        contact_f = Path(pdb_f.replace(".pdb", ".contacts"))
        s = read_pdb(pdb_f)
        clist = get_intermolecular_contacts(s, self.contact_cutoff)
        write_contacts(clist, contact_f)
        if os.path.isfile(contact_f):
            return contact_f
        else:
            return ""

    def cluster(self):
        """Use FCC to cluster structures."""
        ga_log.info(f"[FCC] Calculating contacts, cutoff={self.contact_cutoff}A")
        input_structure_l = self.structure_list[: self.clust_n]
        for pdb in input_structure_l:
            self.calc_contact(pdb)

        contact_file_l = []
        for pdb in input_structure_l:
            contact_f = Path(pdb.replace(".pdb", ".contacts"))
            if contact_f.exists():
                contact_file_l.append(contact_f)

        if not contact_file_l:
            ga_log.warning("No contacts were calculated")

        # Calculate matrix
        ga_log.info("[FCC] Calculating matrix")

        clist = [list(read_contacts(f)) for f in contact_file_l]
        selector1 = selector2 = BY_RESIDUE
        unique = True
        clist_hashed = [
            hash_many(c, unique=unique, selector1=selector1, selector2=selector2)
            for c in clist
        ]
        try:
            idxs, sims = build_matrix(clist_hashed, metric=FCC_METRIC)
        except ZeroDivisionError:
            # Could not built the matrix
            # TODO: check FCC to find out why this happens
            ga_log.warning("[FCC] Could not build the matrix, skipping clustering.")
            return {}
        else:
            if idxs.size == 0:
                ga_log.warning("[FCC] Could not build the matrix, skipping clustering.")
                return {}

        # TODO: Implement a way to re-read the matrix in case its already there
        # idxs, sims = read_matrix(fcc_matrix_f)

        fcc_matrix_f = Path(f"{self.analysis_path}/fcc.matrix")
        write_matrix(idxs, idxs, fcc_matrix_f)

        # cluster
        ga_log.info(
            f"[FCC] Clustering with cutoff={self.clust_cutoff} and minsize={self.clust_min_size}"
        )
        labels = DTB_ALGO(
            idxs, sims, eps=self.clust_cutoff, minsize=self.clust_min_size
        )

        if labels:
            ga_log.info(f"[FCC] {len(labels)} clusters identified")
            ga_log.info("[FCC] Saving cluster information to analysis/cluster.out")
            cluster_out = f"{self.analysis_path}/cluster.out"
            with open(cluster_out, "w") as fh:
                clusters = collections.defaultdict(list)
                for k, v in labels.items():
                    clusters[v].append(k)

                # Sort clusters by size
                sorted_keys = sorted(
                    clusters, key=lambda k: len(clusters.get(k)), reverse=True
                )

                for i, k in enumerate(sorted_keys, start=1):
                    elements = map(str, sorted(clusters[k]))
                    c_str = f'Cluster {i} -> {" ".join(elements)}{os.linesep}'
                    fh.write(c_str)

            # read it again!
            with open(cluster_out, "r") as fh:
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
            ga_log.info("[FCC] No clusters identified")

        return self.cluster_dic

    def evaluate_irmsd(self):
        """Use PROFIT to calculate the interface rmsd."""
        if self.native:
            try:
                ga_log.info("[PROFIT] Evaluating I-RMSD")
                profit = Profit(
                    ref=self.native, mobi=self.structure_list, nproc=self.nproc
                )
                self.irmsd_dic = profit.calc_irmsd()
            except Exception as e:
                ga_log.warning(e)
                ga_log.warning("[PROFIT] Skipping evaluation")
        else:
            ga_log.info(
                "[PROFIT] Native structure not defined, skipping I-RMSD evaluation."
            )

    def output(self):
        """Generate a script friendly output table."""
        output_f = f"{self.analysis_path}/gdock.dat"
        ga_log.info(f"Saving output file to {Path(output_f).name}")

        sep = "\t"
        header = "gen" + sep
        header += "ind" + sep
        header += "ranking" + sep
        header += "satisfaction" + sep
        header += "energy" + sep
        header += "irmsd" + sep
        header += "cluster_id" + sep
        header += "internal_cluster_ranking" + sep
        header += "structure_id" + os.linesep

        with open(output_f, "w") as fh:
            fh.write(header)
            for gen in self.result_dic:
                for ind in self.result_dic[gen]:

                    satisfaction = self.result_dic[gen][ind]["satisfaction"]
                    energy = self.result_dic[gen][ind]["energy"]

                    ranking = self.result_dic[gen][ind]["ranking"]

                    generation_str = str(gen).rjust(4, "0")
                    individual_str = str(ind).rjust(4, "0")

                    structure = self.result_dic[gen][ind]["structure"]
                    if structure:
                        structure_id = pathlib.Path(structure).stem
                        try:
                            irmsd = self.irmsd_dic[structure_id]
                        except KeyError:
                            irmsd = float("nan")
                    else:
                        irmsd = float("nan")

                    # Add the clustering information
                    # TODO: Improve this part, its very messy.
                    model = f"{generation_str}_{individual_str}"
                    internal_ranking = float("nan")
                    cluster_id = float("nan")
                    if self.cluster_dic:
                        for cluster_id in self.cluster_dic:
                            subcluster = self.cluster_dic[cluster_id]
                            if model in subcluster:
                                # the center has already been removed
                                internal_ranking = subcluster.index(model) + 1
                                break
                    if math.isnan(internal_ranking):
                        cluster_id = float("nan")

                    output_str = f"{gen}" + sep
                    output_str += f"{ind}" + sep
                    output_str += f"{ranking}" + sep
                    output_str += f"{satisfaction:.3f}" + sep
                    output_str += f"{energy:.3f}" + sep
                    output_str += f"{irmsd:.2f}" + sep
                    output_str += f"{cluster_id}" + sep
                    output_str += f"{internal_ranking}" + sep
                    output_str += f"{model}" + os.linesep
                    ga_log.debug(output_str)

                    fh.write(output_str)

        fh.close()

        ga_log.info("Generating plots")
        try:
            self.generate_plots(output_f)
        except Exception as e:
            ga_log.debug(e)
            ga_log.warning("Could not generate plots")

    @staticmethod
    def generate_plots(data_f):
        """Generate relevant plots for the data."""
        import warnings

        warnings.filterwarnings("ignore")
        sns.set_context("talk")
        df = pd.read_csv(data_f, sep="\t")
        path = Path(data_f).parent
        kde_outf = Path(path, "kde")
        ridge_outf = Path(path, "ridge")

        # Generate ridgeplot
        sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
        sub_df = df[["gen", "fitness"]]
        quantile_75 = sub_df.quantile([0.75])["fitness"].values[0]
        sub_df = sub_df[sub_df["fitness"] <= quantile_75]

        # jitter the dataframe since there are values that are equal
        #  keep this in mind when presenting the results!
        new_df_list = []
        for e in sub_df.values.tolist():
            gen, value = e
            value += random.uniform(0, 1)
            new_df_list.append((gen, value))

        sub_df = pd.DataFrame(new_df_list, columns=["gen", "fitness"])
        sub_df["gen"] = pd.to_numeric(sub_df["gen"], downcast="integer")

        n_gen = len(set(sub_df["gen"]))
        pal = sns.cubehelix_palette(n_gen)
        g = sns.FacetGrid(
            sub_df, row="gen", hue="gen", aspect=15, height=0.5, palette=pal
        )
        g.map(
            sns.kdeplot,
            "fitness",
            bw_adjust=0.5,
            clip_on=False,
            fill=True,
            alpha=1,
            linewidth=1.5,
        )

        g.map(sns.kdeplot, "fitness", clip_on=False, color="w", lw=2, bw_adjust=0.5)
        g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

        def label(x, color, label):
            ax = plt.gca()
            ax.text(
                0,
                0.2,
                label,
                fontweight="bold",
                color=color,
                ha="left",
                va="center",
                transform=ax.transAxes,
            )

        g.map(label, "fitness")
        g.figure.subplots_adjust(hspace=-0.25)
        g.set_titles("")
        g.set(yticks=[], ylabel="")
        g.despine(bottom=True, left=True)
        g.fig.subplots_adjust(top=0.9)
        g.fig.suptitle("Convergence of fitness over generations")
        g.savefig(ridge_outf.with_suffix(".png"), dpi=300)
        g.savefig(ridge_outf.with_suffix(".svg"))

        # Generate KDE plots
        #  if there is an i-rmsd col = if not all are null
        if not df["irmsd"].isnull().all():
            sns.set_style("white")
            g = sns.jointplot(
                data=df, x="irmsd", y="fitness", kind="kde", color="r", n_levels=5
            )
            g.fig.subplots_adjust(top=0.9)
            g.fig.suptitle("Fitness x Energy distributions")
            g.fig.tight_layout()
            g.savefig(kde_outf.with_suffix(".png"), dpi=300)
            g.savefig(kde_outf.with_suffix(".svg"))
