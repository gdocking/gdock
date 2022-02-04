# A wrapper script to generate ensembles of BM5 using Prody
from pathlib import Path
import os
from prody import parsePDB, ClustENM


def run_clustenm(pdb_f):
    pdb = parsePDB(pdb_f)
    clustenm = ClustENM()
    clustenm.setAtoms(pdb)
    clustenm.writePDBFixed()
    clustenm.run(n_modes=3, n_gens=5,
                 maxclust=tuple(range(20, 120, 20)),
                 sim=True, solvent='imp',
                 t_steps_i=0, t_steps_g=0,
                 platform='CUDA')
    clustenm.writeParameters()
    clustenm.writePDB()


bm5_path = Path('/trinity/login/rodrigo/repos/BM5-clean/HADDOCK-ready')
dirs_to_skip = ['scripts', 'ana_scripts']

for path in bm5_path.glob('*'):
    print(f'> Running for {path}')
    if not path.is_dir():
        continue
    os.chdir(path)
    receptor = path.name + '_r_u.pdb'
    ligand = path.name + '_l_u.pdb'

    receptor_ens = Path(path.name + '_r_u_clustenm.pdb')
    ligand_ens = Path(path.name + '_l_u_clustenm.pdb')

    if not receptor_ens.exists():
        try:
            run_clustenm(receptor)
            print('>> Receptor Done')
        except Exception as e:
            print(e)
            print('>> Receptor Failed')
    if not ligand_ens.exists():
        try:
            run_clustenm(ligand)
            print('>> Ligand Done')
        except Exception as e:
            print(e)
            print('>> Ligand Failed')

    print(f'>> {path} Done')