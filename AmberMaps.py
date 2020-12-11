__doc__='''
This module include mapping required by the workflow.
------------------------------------------------------------
(MAP) Resi_map
Map residue key with the 3 letter name.
{resi_key:resi_name_3,...}
------------------------------------------------------------
(MAP) Resi_map2
Map 3 letter name with residue key.
{resi_name_3:resi_key,...}
------------------------------------------------------------
(List) Resi_list
Map residue key with index number
[resi_key,...]
------------------------------------------------------------
(???) Atom_topology


------------------------------------------------------------
'''

Resi_map={'R':'ARG',
          'H':'HIS',
          'K':'LYS',
          'D':'ASP',
          'E':'GLU',
          'S':'SER',
          'T':'THR',
          'N':'ASN',
          'Q':'GLN',
          'C':'CYS',
          'U':'SEC',
          'G':'GLY',
          'P':'PRO',
          'A':'ALA',
          'V':'VAL',
          'I':'ILE',
          'L':'LEU',
          'M':'MET',
          'F':'PHE',
          'Y':'TYR',
          'W':'TRP'}

Resi_map2={'ARG':'R',
           'HIS':'H',
           'HIE':'H',
           'HIP':'H',
           'HID':'H',
           'LYS':'K',
           'ASP':'D',
           'GLU':'E',
           'SER':'S',
           'THR':'T',
           'ASN':'N',
           'GLN':'Q',
           'CYS':'C',
           'SEC':'U',
           'GLY':'G',
           'PRO':'P',
           'ALA':'A',
           'VAL':'V',
           'ILE':'I',
           'LEU':'L',
           'MET':'M',
           'PHE':'F',
           'TYR':'Y',
           'TRP':'W'}

Resi_list=['R','H','K','D','E','S','T','N','Q','C','U','G','P','A','V','I','L','M','F','Y','W']

TIP3P_map=[ "O" , "H" , "OW", "HW","HOH","WAT",
            "F" , 
	        "Cl", "CL",
	        "Br", "BR",
	        "I" ,
	        "Li", "LI",
	        "Na", "NA",
            "K" ,
	        "Rb", "RB",
	        "Cs", "CS",
	        "Mg", "MG",
            "Tl", "TL",
            "Cu", "CU",
            "Ag", "AG",
            "Be", "BE",
            "Cu", "CU",
            "Ni", "NI",
            "Pt", "PT",
            "Zn", "ZN",
            "Co", "CO",
            "Pd", "PD",
            "Ag", "AG",
            "Cr", "CR",
            "Fe", "FE",
            "Mg", "MG",
            "V" ,
            "Mn", "MN",
            "Hg", "HG",
            "Cd", "CD",
            "Yb", "YB",
            "Ca", "CA",
            "Sn", "SN",
            "Pb", "PB",
            "Eu", "EU",
            "Sr", "SR",
            "Sm", "SM",
            "Ba", "BA",
            "Ra", "RA",
            "Al", "AL",
            "Fe", "FE",
            "Cr", "CR",
            "In", "IN",
            "Tl", "TL",
            "Y" ,
            "La", "LA",
            "Ce", "CE",
            "Pr", "PR",
            "Nd", "ND",
            "Sm", "SM",
            "Eu", "EU",
            "Gd", "GD",
            "Tb", "TB",
            "Dy", "DY",
            "Er", "ER",
            "Tm", "TM",
            "Lu", "LU",
            "Hf", "HF",
            "Zr", "ZR",
            "Ce", "CE",
            "U" ,
            "Pu", "PU",
            "Th", "TH"  ]
