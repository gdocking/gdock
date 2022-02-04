import logging
from pdbtools.pdb_splitmodel import split_model
import tempfile
from pathlib import Path
import os
import glob
from gdock.modules.error import PDBError

ga_log = logging.getLogger("ga_log")


class PDB:
    def __init__(self):
        """Initialize PDB class."""
        self.coords = {}
        self.raw_pdb = {}

    def load(self, pdb_f):
        """Load a PDB file into a dictionary."""
        ga_log.debug(f"Loading {pdb_f}")

        chain_list = self._identify_chains(pdb_f)

        if len(chain_list) > 1:
            raise PDBError(
                f"{pdb_f} contains multiple chains {chain_list}, we are expecteding one chain per PDB."
            )

        target_chain = chain_list[0]
        if target_chain in self.raw_pdb:
            raise PDBError(
                f"{pdb_f} contains chain {target_chain} which has already been loaded."
            )

        self.raw_pdb[target_chain] = {}
        self.coords[target_chain] = {}

        current_dir = Path.cwd()
        with tempfile.TemporaryDirectory() as tmpdirname:
            # The input might be an ensemble, use a temporary directory
            #  to keep things organized
            os.chdir(tmpdirname)

            with open(pdb_f, "r") as fh:
                split_model(fh)

            pdb_list = glob.glob(f"{tmpdirname}/*pdb")
            if not pdb_list:
                # There were no models, this means the initial pdb
                #  is not an ensemble
                pdb_list = [pdb_f]

            for i, pdb in enumerate(pdb_list):

                self.raw_pdb[target_chain][i] = []
                self.coords[target_chain][i] = []

                with open(pdb, "r") as fh:
                    for line in fh.readlines():
                        if line.startswith("ATOM"):
                            chain = line[21]
                            x = float(line[31:38])
                            y = float(line[39:46])
                            z = float(line[47:54])
                            self.coords[chain][i].append((x, y, z))
                            self.raw_pdb[chain][i].append(line)
            fh.close()
        os.chdir(current_dir)

    @staticmethod
    def _identify_chains(pdb_f):
        """Identify which chains are present in the PDB."""
        chain_l = []
        with open(pdb_f, "r") as fh:
            for line in fh.readlines():
                if line.startswith("ATOM"):
                    chain = line[21]
                    chain_l.append(chain)
        return list(set(chain_l))


class Restraint:
    def __init__(self, raw_pdb):
        """Initialize Restraint class."""
        self.raw_pdb = raw_pdb
        self.coords = {}
        self.resnums = {}

    def load(self, restraint, identifier):
        """Identity the coordinates of the restraints residue numbers."""
        ga_log.debug(f"Loading restraints {identifier}, {restraint}")
        self.resnums[identifier] = restraint
        if identifier not in self.coords:
            self.coords[identifier] = {}

        for model in self.raw_pdb[identifier]:
            # for chain in self.raw_pdb:
            # if chain == identifier:
            self.coords[identifier][model] = []
            for line in self.raw_pdb[identifier][model]:
                resnum = int(line[22:26])
                if resnum in restraint:
                    x = float(line[31:38])
                    y = float(line[39:46])
                    z = float(line[47:54])
                    self.coords[identifier][model].append((x, y, z))
