from subprocess import run, CalledProcessError

target_dir = './Input_files_SI'
tags= []
with open('Mutants.csv') as f:
    for line in f:
        tags.append(line.strip())
count = 0
for i in tags:
    count += 1
    search_name_inp = 'FAcD-FA-ASP_rmW_'+i+'_rmH_aH_min_rmW.inpcrd'
    search_name_prm = 'FAcD-FA-ASP_rmW_'+i+'_rmH_aH_min_rmW.prmtop'
    inp_path = run('find . -name '+search_name_inp, check=True, text=True, shell=True, capture_output=True).stdout.strip()
    prm_path = run('find . -name '+search_name_prm, check=True, text=True, shell=True, capture_output=True).stdout.strip()
    cmd = ' '.join(["cp",inp_path,prm_path,target_dir])
    print(count)
    try:
        run(cmd, check=True, text=True, shell=True, capture_output=True)
    except CalledProcessError:
        print(i)
    

    
