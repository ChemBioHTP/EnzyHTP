from math import frexp

from core.clusters.accre import Accre
from Class_PDB import PDB
from Class_Conf import Config
from helper import write_data

# settings
Config.n_cores = 4
Config.max_core = 2000
Config.Amber.conf_equi['nstlim'] = 500000
Config.Amber.conf_prod['nstlim'] = 55000000
Config.Amber.conf_prod['ntwx'] = '50000'
Config.debug = 1
data_output_path = './Mutation.dat'

def main():
    # Mutation
    start_variant = "KE-07_aH_IA7D_KA146E_GA202R_NA224D.pdb"
    selected_mutants = [["V12M", "E146K"], "L47I"]
    for i, mutaflags in enumerate(selected_mutants):
        pdb_obj = PDB(start_variant, wk_dir=f"./mutation_{i}")
        pdb_obj.Add_MutaFlag(mutaflags)
        pdb_obj.PDB2PDBwLeap()
        # use minimization to relax each mutated PDB
        pdb_obj.PDB2FF(local_lig=0, ifsavepdb=1)
        pdb_obj.conf_min['DISANG']=f'{pdb_obj.cache_path}/PDBMin/0.rs'
        pdb_obj.conf_min['nmropt_rest'] = '1'
        pdb_obj.get_stru()
        a1=str(pdb_obj.stru.ligands[0].CAE.id)
        a2=str(pdb_obj.stru.ligands[0].H2.id)
        a3=str(pdb_obj.stru.chains[0].i101.OE2.id)
        pdb_obj.conf_min['rs_constraints']=[
                                {
                                'iat':[a2,a3],'r1': '0.50','r2': '0.96','r3': '2.50','r4': '3.50',
                                'rk2':'0.0','rk3':'100.0','ir6':'1','ialtd':'0'
                                },
                                {
                                'iat':[a1,a2,a3],'r1': '120.0','r2': '150.0','r3': '210.0','r4': '250.0',
                                'rk2':'100.0','rk3':'100.0','ir6':'1','ialtd':'0'
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
        pdb_obj.PDB2FF(local_lig=0, ifsavepdb=1) # need to update stru be for getting
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
                        period=30,
                        res_setting={'account':'csb_gpu_acc'} )

        # QM
        pdb_obj.nc2mdcrd(start=101,step=10)
        atom_mask = ':101,254'
        g_route = '# pbe1pbe/def2SVP nosymm'
        pdb_obj.PDB2QMCluster(  atom_mask, 
                                g_route=g_route,
                                ifchk=1,
                                if_cluster_job=1, 
                                cluster=Accre(), 
                                job_array_size=20,
                                period=30,
                                res_setting={'account':'yang_lab'} )
        pdb_obj.get_fchk(keep_chk=0)

        # --- Analysis ---
        # targeting C-I bond
        a1 = int(pdb_obj.stru.ligands[0].CAE)
        a2 = int(pdb_obj.stru.ligands[0].H2)
        a1qm = pdb_obj.qm_cluster_map[str(a1)]
        a2qm = pdb_obj.qm_cluster_map[str(a2)]
        # Field Strength (MM)
        E_atom_mask = ':1-100,102-253'
        Es = pdb_obj.get_field_strength(E_atom_mask, a1=a1 ,a2=a2 ,bond_p1='center') 
        # Bond Dipole Moment (QM)
        Dipoles = PDB.get_bond_dipole(pdb_obj.qm_cluster_fchk, a1qm, a2qm)

        # # SASA
        # sasa_mask = ""
        # sasa = pdb_obj.get_sasa(sasa_mask)
        # # RMSD
        # rmsd_mask = ""
        # rmsd = pdb_obj.get_rmsd(rmsd_mask)

        write_data(pdb_obj.MutaFlags, {'E': Es, 'Bond Dipole': Dipoles, 'traj': pdb_obj.mdcrd},  data_output_path)


if __name__ == "__main__":
    main()
