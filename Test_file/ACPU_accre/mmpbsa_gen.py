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

mutants = XXX
wt_pdb = "YYY"

def main():
    # Mutation
    for mut in mutants:
        pdb_obj = PDB(wt_pdb, wk_dir=f"./mutation_{'_'.join(mut)}")
        pdb_obj.Add_MutaFlag(mut)
        print(f"Working on mutation: {pdb_obj.MutaFlags}")
        pdb_obj.PDB2PDBwLeap()
        # use minimization to relax each mutated PDB
        pdb_obj.PDB2FF(local_lig=0, 
                    ifsavepdb=1, 
                    lig_net_charge_mapper={
                            ZZZ
                    }) # need to update stru be for getting
        pdb_obj.conf_min['DISANG']=f'{pdb_obj.cache_path}/PDBMin/0.rs'
        pdb_obj.conf_min['nmropt_rest'] = '1'
        pdb_obj.get_stru()
        a1=str(pdb_obj.stru.ligands[1].QQQ.id)
        a2=str(pdb_obj.stru.ligands[0].C64.id)
        pdb_obj.conf_min['rs_constraints']=[
                        {
                        'iat':[a1,a2],'r1': '0.50','r2': '2.3','r3': '2.5','r4': '2.7',
                        'rk2':'0.0','rk3':'100.0','ir6':'1','ialtd':'0'
                        },
                        ]
        pdb_obj.PDBMin( engine='Amber_CPU', 
                        if_cluster_job=1,
                        cluster=Accre(),
                        period=10,
                        res_setting={'node_cores': '24',
                                    'mem_per_core' : '3G',
                                    'account':'yang_lab'} )
        # protonation modification
        pdb_obj.rm_wat()
        pdb_obj.rm_allH()
        pdb_obj.get_protonation(if_prt_ligand=0)

        # MD sampling
        disang_path = f'{pdb_obj.dir}/MD/0.rs'
        pdb_obj.conf_min['nmropt_rest'] = '1'
        pdb_obj.conf_min['DISANG']=disang_path
        pdb_obj.conf_heat['nmropt_rest'] = '1'
        pdb_obj.conf_heat['DISANG']=disang_path
        pdb_obj.conf_equi['nmropt_rest'] = '1'
        pdb_obj.conf_equi['DISANG']=disang_path
        pdb_obj.conf_prod['nmropt_rest'] = '1'
        pdb_obj.conf_prod['DISANG']=disang_path
        pdb_obj.PDB2FF(local_lig=0, 
                    ifsavepdb=1, 
                    lig_net_charge_mapper={
                            ZZZ
                    }) # need to update stru be for getting
        pdb_obj.get_stru()
        a1=str(pdb_obj.stru.ligands[1].QQQ.id)
        a2=str(pdb_obj.stru.ligands[0].C64.id)
        constrain_dict=[
                        {
                        'iat':[a1,a2],'r1': '0.50','r2': '2.3','r3': '2.5','r4': '2.7',
                        'rk2':'0.0','rk3':'100.0','ir6':'1','ialtd':'0'
                        },
                        ]
        pdb_obj.conf_min['rs_constraints']=constrain_dict
        pdb_obj.conf_heat['rs_constraints']=constrain_dict
        pdb_obj.conf_equi['rs_constraints']=constrain_dict
        pdb_obj.conf_prod['rs_constraints']=constrain_dict

        pdb_obj.PDBMD(  engine='Amber_GPU', 
                        equi_cpu=0, 
                        if_cluster_job=1,
                        cluster=Accre(),
                        period=30,
                        res_setting={'account':'csb_gpu_acc'} )
        pdb_obj.nc2mdcrd(start=101,step=10)

        # --- Analysis ---
        # MMPBSA
        ligand_mask = ":902"
        mmpbsa_result_dict = pdb_obj.get_mmpbsa_binding(
            ligand_mask,
            cluster=Accre(),
            res_setting = {'account':'yang_lab'})

        write_data(
            pdb_obj.MutaFlags, 
            {'gb_binding': mmpbsa_result_dict['gb']['mean']['DELTA TOTAL'], 'pb_binding': mmpbsa_result_dict['pb']['mean']['DELTA TOTAL']},
            data_output_path)


if __name__ == "__main__":
    main()
