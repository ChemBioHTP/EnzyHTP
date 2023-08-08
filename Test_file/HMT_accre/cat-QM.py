import datetime
from glob import glob
from Class_PDB import *
from Class_Conf import *
from Class_ONIOM_Frame import *
from helper import write_data, line_feed

# settings
Config.n_cores = 12
Config.max_core = 2000
#Config.PC_cmd = 'srun' # srun does not work like this use mpi instead
Config.debug = 1
data_output_path = './Mutation.dat'

def main():

	pdb_obj = PDB(glob('*_ff.pdb')[0])
	pdb_obj.Add_MutaFlag('XXX')
	pdb_obj.prmtop_path = glob('*prmtop')[0]
	pdb_obj.prepi_path = {'SAH':'../ligands/ligand_SAH.prepin', 'MET':'../ligands/ligand_MET.prepin'}
	pdb_obj.mdcrd = glob('./MD/*mdcrd')[0]
	# --- QM cluster ---
	atom_mask = ':217,218'
	pdb_obj.get_stru()
	sele_lines, pdb_obj.qm_cluster_map = pdb_obj.stru.get_sele_list(atom_mask, fix_end='H', prepi_path=pdb_obj.prepi_path)
	pdb_obj.qm_cluster_fchk = glob('./QMCluster/*fchk')
	pdb_obj.qm_cluster_fchk.sort(key=lambda file_name: int(file_name.split('.')[-2].split('_')[-1]))

	# pdb_obj.get_fchk(keep_chk=0)

	# --- Analysis ---
	# targeting C-I bond
	a1 = int(pdb_obj.stru.ligands[1].C1)
	a2 = int(pdb_obj.stru.ligands[1].I1)
	a1qm = pdb_obj.qm_cluster_map[str(a1)]
	a2qm = pdb_obj.qm_cluster_map[str(a2)]
	# Field Strength (MM)
	E_atom_mask = ':1-216'
	Es = pdb_obj.get_field_strength(E_atom_mask, a1=a1 ,a2=a2 ,bond_p1='center') 
	# Bond Dipole Moment (QM)
	print(a1qm, a2qm)
	print(pdb_obj.qm_cluster_fchk)
	Dipoles = PDB.get_bond_dipole(pdb_obj.qm_cluster_fchk, a1qm, a2qm)
	# write to csv or plot
	write_data(pdb_obj.MutaFlags, {'E': Es, 'Bond Dipole': Dipoles},  data_output_path)

if __name__ == "__main__":
	main()
