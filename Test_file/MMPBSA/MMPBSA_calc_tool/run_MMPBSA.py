import os
out_dir = '/scratch/shaoq1/CHEM5420-PPIranking/MMPBSA/'

def main():
    with open('MMPBSA_input.list') as f:
        data_paths = eval(f.read().strip())
        for inp_list in data_paths:
            print('Working on: '+inp_list[0])
            out_path = out_dir+inp_list[0]+'.dat'
            nc = inp_list[1]
            prmtop = inp_list[2]
            dl_prmtop = inp_list[3]
            dr_prmtop = inp_list[4]
            dc_prmtop = inp_list[5]
            cmd = 'srun python $AMBERHOME/bin/MMPBSA.py.MPI -O -i /scratch/shaoq1/CHEM5420-PPIranking/MMPBSA/MMPBSA.in -o '+out_path+' -sp '+ prmtop +' -cp '+ dc_prmtop +' -rp '+ dr_prmtop +' -lp '+ dl_prmtop + ' -y '+ nc
            print('running:  '+cmd)
            os.system(cmd)
            os.system('rm reference.frc')


if __name__ == "__main__":
    main()
