f1 = open('lasso_extension/scaffolds/8mr_lassoPeptide_asx_8A.pdb')
f2 = open('lasso_extension/output_structures/8mr_asx_8A.pdb')

lines1 = f1.readlines()[1:]
lines2 = f2.readlines()[2:]
f1.close()
f2.close()
for i in range(len(lines1)):
    l1 = lines1[i].split()[6:9]
    l2 = lines2[i].split()[5:8]
    if l1 != l2:
        print(lines1[i])