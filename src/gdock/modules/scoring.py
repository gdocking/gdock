import configparser
import logging
import multiprocessing
import os
import pathlib
import shlex
import subprocess  # nosec

from gdock.modules.files import get_full_path

ga_log = logging.getLogger("ga_log")

etc_folder = get_full_path("etc")
ini = configparser.ConfigParser(os.environ)
ini.read(os.path.join(etc_folder, "gdock.ini"), encoding="utf-8")
dcomplex_exe = ini.get("third_party", "dcomplex_exe")


class Scoring:
    def __init__(self, data_dic, params):
        self.dcomplex_exec = dcomplex_exe
        self.data_dic = data_dic
        self.nproc = params["main"]["number_of_processors"]
        self.structure_list = []
        self.scoring_dic = {}
        self.ranked = {}

    def score(self):
        """Score the structures."""
        ga_log.info("Scoring structures")
        self.scoring_dic = {}
        for gen in self.data_dic:
            for ind in self.data_dic[gen]:
                pdb_f = self.data_dic[gen][ind]["structure"]
                if pdb_f:
                    fitness = self.data_dic[gen][ind]["fitness"]
                    self.structure_list.append(pdb_f)

                    pdb_name = pathlib.Path(pdb_f).stem

                    satisfaction = fitness[0]
                    # other_fitness = fitness[1]
                    self.scoring_dic[pdb_name] = [satisfaction]

        # return ordered_list_of_ranked_structures
        energy_structure_dic = self.calculate_energy()
        for pdb_name in energy_structure_dic:
            energy_v = energy_structure_dic[pdb_name]
            self.scoring_dic[pdb_name].append(energy_v)

        # magic time, do the scoring!
        # gdock score?
        ranked_list = []
        for pdb_id in self.scoring_dic:
            satisfaction = self.scoring_dic[pdb_id][0]
            energy = self.scoring_dic[pdb_id][1]

            try:
                # this weights were obtained via optimize_score.py,
                #  but did not improve the success rate, keep it here anyway
                # gdock_score = (energy * -0.38) / (satisfaction * -0.49)
                gdock_score = energy / satisfaction

            except ZeroDivisionError:
                gdock_score = 0.0

            if gdock_score != 0.0:
                ranked_list.append((pdb_id, gdock_score))

        counter = 1
        for pdb_id, score in sorted(ranked_list, key=lambda tup: tup[1]):
            self.ranked[pdb_id] = {"rank": counter, "score": score}
            counter += 1

        # add the ranking to the data structure
        for gen in self.data_dic:
            for ind in self.data_dic[gen]:
                pdb_f = self.data_dic[gen][ind]["structure"]
                ranking = float("nan")
                score = float("nan")
                energy = float("nan")

                if pdb_f:

                    pdb_id = pathlib.Path(pdb_f).stem
                    energy = self.scoring_dic[pdb_id][1]

                    try:
                        ranking = self.ranked[pdb_id]["rank"]
                        score = self.ranked[pdb_id]["score"]
                    except KeyError:
                        # no ranking/score for these, check above what is the
                        #  criteria for it to be included in self.ranked; it
                        #  could be that it has a structure but
                        #  its satisfaction is 0
                        pass

                self.data_dic[gen][ind]["ranking"] = ranking
                self.data_dic[gen][ind]["score"] = score
                self.data_dic[gen][ind]["energy"] = energy

        return self.data_dic

    def calculate_energy(self):
        """Obtain energies for the internal list of structures."""
        energy_dic = {}
        pool = multiprocessing.Pool(processes=self.nproc)
        results = []
        for structure in self.structure_list:
            results.append(
                (
                    structure,
                    pool.apply_async(
                        self.run_dcomplex, args=(self.dcomplex_exec, structure)
                    ),
                )
            )

        for p in results:
            pdb, irmsd = p[0], p[1].get()
            pdb_identifier = pathlib.Path(pdb).stem
            energy_dic[pdb_identifier] = irmsd

        pool.close()
        pool.join()

        return energy_dic

    @staticmethod
    def run_dcomplex(dcomplex_exe, pdb_f):
        """Use DCOMPLEX to calculate the PDB's energy."""
        cmd = f"{dcomplex_exe} {pdb_f} A B"
        ga_log.debug(f"cmd is: {cmd}")
        out = subprocess.check_output(shlex.split(cmd), shell=False)  # nosec
        result = out.decode("utf-8").split("\n")
        energy = float(result[-2].split()[1])
        return energy
