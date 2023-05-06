from glob import glob
from core.clusters.accre import Accre
from Class_PDB import PDB
from Class_Conf import Config
from helper import write_data

# settings
Config.n_cores = 1
Config.max_core = 2000
Config.Amber.conf_equi['nstlim'] = 500000
Config.Amber.conf_prod['nstlim'] = 55000000
Config.Amber.conf_prod['ntwx'] = '50000'
Config.debug = 1
data_output_path = './Mutation.dat'

mut = XXX
wt_pdb = "puo_YYY.pdb"

def main():
    pdb_obj = PDB(wt_pdb, wk_dir=f"./mutation_{'_'.join(mut)}")
    pdb_obj.Add_MutaFlag(mut)
    print(f"Working on mutation: {pdb_obj.MutaFlags}")

    pdb_obj.nc = f"{pdb_obj.dir}/MD/prod.nc"
    pdb_obj.prmtop_path = glob(f"{pdb_obj.dir}/*_ff_min_rmW_rmH_aH.prmtop")[0]
    pdb_obj.path = glob(f"{pdb_obj.dir}/*_ff_min_rmW_rmH_aH_ff.pdb")[0]

    pdb_obj.nc2mdcrd(start=101,step=10)
    # --- Analysis ---
    # MMPBSA
    ligand_mask = ":902"
    mmpbsa_result_dict = pdb_obj.get_mmpbsa_binding(
        ligand_mask,
        cluster=Accre(),
        res_setting = {'account':'yang_lab_csb'})

    write_data(
        pdb_obj.MutaFlags, 
        {'gb_binding': mmpbsa_result_dict['gb']['mean']['DELTA TOTAL'], 'pb_binding': mmpbsa_result_dict['pb']['mean']['DELTA TOTAL']},
        data_output_path)


if __name__ == "__main__":
    main()
