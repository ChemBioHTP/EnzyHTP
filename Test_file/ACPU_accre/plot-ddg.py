
#Set project path here
proj_path='./'
csv_file = f'{proj_path}PuO_1st_ddg.csv'
resi_idx_mapper = {
    "TynA":{
        "A": -5,
        "B": 715
        },
    "PuO":{
        "A": -1,
        "B": 449 
    }
}

def extract_data_ddg():
    """
    extract ddg data from ddg_result.dat
    """
    result = []
    with open("ddg_results.dat") as f:
        ddg_mapper = eval(f.read().strip())
    for k, v in ddg_mapper.items():
        new_k = []
        existing_r_idx = []
        for k_p in k:
            orig, r_idx, target = k_p.strip().split()
            if r_idx not in existing_r_idx and str(int(r_idx)+resi_idx_mapper["PuO"]["A"]-resi_idx_mapper["PuO"]["B"]) not in existing_r_idx:
                existing_r_idx.append(r_idx)
                new_k.append("".join((orig, r_idx, target)))
        new_k = " ".join(new_k)
        result.append((new_k, v))

    return result

data = extract_data_ddg()

# Output
#=========
#---csv---
with open(csv_file,'w') as of:
    # title
    of.write("Mutant,cart_ddg\n")
    for m_name, m_data in data:
        of.write(f"{m_name},{m_data}\n")

