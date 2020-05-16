# GA needs to solve the Quaternion (grid position implemented later)
from functions import calc_clash
from modules.ga import GeneticAlgorithm
from modules.structure import PDB

if __name__ == '__main__':
    # 1. Load structure
    pdb = PDB()
    pdb.load('target-unbound.pdb')
    # 2. Position
    pass
    # DEV: randomize orientation of chain B
    pdb.randomize_rotation('B')
    pdb.output('rotated.pdb')
    calc_clash('rotated.pdb')

    ## Proof-of-concept
    # from a unbound complex, randomize the ligand initial orientation and use GA to find it again, use clash as fitness
    ga = GeneticAlgorithm('rotated.pdb', population_size=100, number_of_generations=30, target_chain='B', nproc = 4)
    toolbox = ga.setup()
    result_dic = ga.run(toolbox)
