from HTP_Traj_calc import *
from Class_Conf import *


def main():
    # Config.n_cores=24  
    #Config.PC_cmd='srun'  
    #Config.Amber.MMPBSA.MMPBSA_EXE='/scratch/shaoq1/CHEM5420-PPIranking/MMPBSA.py.MPI'
    igb=5
    frag_str = 'c. A+B:c. C'

    traj_list = Traj_calc.From_keyword_find('.', keyword='1T83*')
    input_files = Traj_calc.calc_MMPBSA(traj_list, frag_str, igb, use_parmed=0, prepare_only=1)
    with open('MMPBSA_input.list','w') as of:
        of.write(repr(input_files))


if __name__ == "__main__":
    main()
