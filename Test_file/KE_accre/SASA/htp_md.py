from core.clusters.accre import Accre
from Class_PDB import PDB
from Class_Conf import Config
from helper import write_data, check_complete_metric_run

from enzy_htp.core.file_system import get_valid_temp_name

# settings
Config.n_cores = 1
Config.max_core = 2000
Config.Amber.conf_equi['nstlim'] = 50000 # note now its 1:11 instead of 1:110
Config.Amber.conf_prod['nstlim'] = 550000
Config.Amber.conf_prod['ntwx'] = '500'
Config.debug = 1
data_output_path = './Mutation.dat'

wt_pdb = "YYY"

def main():
    # Mutation
    for i in range(10):
        pdb_obj = PDB(wt_pdb, wk_dir=f"./rep_{i}")
        # protonation modificationm
        pdb_obj.rm_allH()
        pdb_obj.get_protonation(if_prt_ligand=0)

        # MD sampling
        disang_path = get_valid_temp_name('./0_md.rs')
        pdb_obj.conf_min['nmropt_rest'] = '1'
        pdb_obj.conf_min['DISANG']=disang_path
        pdb_obj.conf_heat['nmropt_rest'] = '1'
        pdb_obj.conf_heat['DISANG']=disang_path
        pdb_obj.conf_equi['nmropt_rest'] = '1'
        pdb_obj.conf_equi['DISANG']=disang_path
        pdb_obj.conf_prod['nmropt_rest'] = '1'
        pdb_obj.conf_prod['DISANG']=disang_path
        pdb_obj.PDB2FF(local_lig=0, 
                    ifsavepdb=1, ) # need to update stru be for getting
        pdb_obj.get_stru()
        a1=str(pdb_obj.stru.ligands[0].CAE.id)
        a2=str(pdb_obj.stru.ligands[0].H2.id)
        a3=str(pdb_obj.stru.chains[0].i101.OE2.id)
        constrain_dict = [{
                            'iat':[a2,a3],'r1': '0.50','r2': '0.96','r3': '2.00','r4': '2.50',
                            'rk2':'0.0','rk3':'100.0','ir6':'1','ialtd':'0'
                            },
                            {
                            'iat':[a1,a2,a3],'r1': '120.0','r2': '150.0','r3': '210.0','r4': '250.0',
                            'rk2':'100.0','rk3':'100.0','ir6':'1','ialtd':'0'
                            },]
        pdb_obj.conf_min['rs_constraints']=constrain_dict
        pdb_obj.conf_heat['rs_constraints']=constrain_dict
        pdb_obj.conf_equi['rs_constraints']=constrain_dict
        pdb_obj.conf_prod['rs_constraints']=constrain_dict

        pdb_obj.PDBMD(  engine='Amber_GPU', 
                        equi_cpu=0, 
                        if_cluster_job=1,
                        cluster=Accre(),
                        period=600,
                        res_setting={'account':'csb_gpu_acc'} )
        pdb_obj.nc2mdcrd(start=101,step=10) # sample per 10 ps

        # --- Analysis ---
        a1 = int(pdb_obj.stru.ligands[0].CAE)
        a2 = int(pdb_obj.stru.ligands[0].H2)
        E_atom_mask = ':1-100,102-253'
        Es = pdb_obj.get_field_strength(E_atom_mask, a1=a1 ,a2=a2 ,bond_p1='center') 

        # SASA
        mask_sasa = ":9,11,48,50,101,128,201,202,222"
        mask_pro = ":1-253"
        mask_sub = ":254"
        sasa_value = PDB.get_sasa_ratio(str(pdb_obj.prmtop_path), str(pdb_obj.mdcrd), 
                                        mask_pro, mask_sasa, mask_sub)
        write_data(
            pdb_obj.MutaFlags, 
            {'E': Es, 'sasa_avg': sasa_value},
            data_output_path)

if __name__ == "__main__":
    main()
