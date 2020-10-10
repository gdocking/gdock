import logging

ga_log = logging.getLogger('ga_log')


class PDB:

    def __init__(self):
        """Initialize PDB class."""
        self.coords = {}
        self.raw_pdb = {}

    def load(self, pdb_f):
        """Load a PDB file into a dictionary."""
        ga_log.debug(f'Loading {pdb_f}')
        with open(pdb_f, 'r') as fh:
            for line in fh.readlines():
                if line.startswith('ATOM'):
                    chain = line[21]
                    if chain not in self.coords:
                        self.coords[chain] = []
                    if chain not in self.raw_pdb:
                        self.raw_pdb[chain] = []
                    x = float(line[31:38])
                    y = float(line[39:46])
                    z = float(line[47:54])
                    self.coords[chain].append((x, y, z))
                    self.raw_pdb[chain].append(line)
        fh.close()


class Restraint:
    def __init__(self, raw_pdb):
        """Initialize Restraint class."""
        self.raw_pdb = raw_pdb
        self.coords = {}

    def load(self, restraint, identifier):
        """Identity the coordinates of the restraints residue numbers."""
        ga_log.debug(f'Loading restraints {identifier}, {restraint}')
        if identifier not in self.coords:
            self.coords[identifier] = []

        for line in self.raw_pdb[identifier]:
            resnum = int(line[22:26])
            if resnum in restraint:
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                self.coords[identifier].append((x, y, z))
