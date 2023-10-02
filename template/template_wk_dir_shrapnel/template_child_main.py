'''The template for the old EnzyHTP main script:
A workflow that calculate properties for mutants based on QM and MM.

Author: QZ Shao <shaoqz@icloud.com>
Date: 2023-8-8'''
import pickle

from core.clusters.accre import Accre
from Class_PDB import PDB
from Class_Conf import Config
from helper import write_data


# Configurations
## resource of the main script
Config.n_cores = 1
Config.max_core = 2000 #MB
## Details of MD
Config.Amber.conf_equi['nstlim'] = 500000 # * 2 fs = 1 ns
Config.Amber.conf_prod['nstlim'] = 55000000 # * 2 fs = 110 ns
Config.Amber.conf_prod['ntwx'] = '50000' # * 2 fs = 0.1 ns

# Input
mutants = XXX
wt_pdb = "YYY"
# Output
data_output_path_dat = './Mutation.dat'


def main():
    for mut in mutants:        
        # Prepare
        pdb_obj = PDB(wt_pdb, wk_dir=f"./mutation_{'_'.join(mut)}")
        pdb_obj.rm_allH()
        pdb_obj.get_protonation(if_prt_ligand=0)

        # Mutation
        pdb_obj.Add_MutaFlag(mut)
        pdb_obj.PDB2PDBwLeap()
        ## use minimization to relax the crude initial mutant structure
        pdb_obj.PDB2FF(local_lig=0, ifsavepdb=1)
        pdb_obj.PDBMin(cycle=20000,
                       engine='Amber_CPU', 
                       if_cluster_job=1,
                       cluster=Accre(),
                       period=180,
                       res_setting={'node_cores': '24',
                                    'mem_per_core' : '3G',
                                    'account':'xxx'} )
        pdb_obj.rm_wat()
        ## protonation perturbed by mutations
        pdb_obj.rm_allH()
        pdb_obj.get_protonation(if_prt_ligand=0)

        # MD sampling
        pdb_obj.PDB2FF(local_lig=0, ifsavepdb=1)
        pdb_obj.PDBMD(engine='Amber_GPU', 
                      if_cluster_job=1,
                      cluster=Accre(),
                      period=600,
                      res_setting={'account':'xxx'} )
        ## sample from traj (.nc file)
        pdb_obj.nc2mdcrd(start=101,step=10)

        # QM Cluster
        atom_mask = ':101,254'
        g_route = '# pbe1pbe/def2SVP nosymm'
        pdb_obj.PDB2QMCluster(  atom_mask, 
                                g_route=g_route,
                                ifchk=1,
                                if_cluster_job=1, 
                                cluster=Accre(), 
                                job_array_size=20,
                                period=120,
                                res_setting={'account':'xxx'} )
        pdb_obj.get_fchk(keep_chk=0)

        # --- Analysis ---
        pdb_obj.get_stru()
        # targeting C-I bond
        a1 = int(pdb_obj.stru.ligands[0].CAE)
        a2 = int(pdb_obj.stru.ligands[0].H2)
        a1qm = pdb_obj.qm_cluster_map[str(a1)]
        a2qm = pdb_obj.qm_cluster_map[str(a2)]
        # Field Strength (MM)
        e_atom_mask = ':1-100,102-253'
        e_list = pdb_obj.get_field_strength(
            e_atom_mask,
            a1=a1, a2=a2, bond_p1='center') 
        # Bond Dipole Moment (QM)
        dipole_list = PDB.get_bond_dipole(pdb_obj.qm_cluster_fchk, a1qm, a2qm)

        # SASA ratio
        mask_sasa = ":9,11,48,50,101,128,201,202,222"
        mask_pro = ":1-253"
        mask_sub = ":254"
        sasa_ratio = PDB.get_sasa_ratio(str(pdb_obj.prmtop_path), str(pdb_obj.mdcrd), 
                                        mask_pro, mask_sasa, mask_sub)

        # write output (readable style)
        write_data(
            pdb_obj.MutaFlags, 
            {
            'field_strength': e_list,
            'bond_dipole': dipole_list,
            'sasa_ratio': sasa_ratio,
            'traj': pdb_obj.mdcrd,
            },
            data_output_path_dat)


if __name__ == "__main__":
    main()
