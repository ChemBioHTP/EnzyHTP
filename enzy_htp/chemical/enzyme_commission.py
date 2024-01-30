"""TODO(CJ)


Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-10-30
"""

#TODO(CJ): Need a lot more documentation here
from rdkit.Chem import AllChem

import numpy as np

from typing import Dict
import pickle
import pandas as pd
import editdistance
import requests
from bs4 import BeautifulSoup

from enzy_htp.core import _LOGGER

BRENDA_MAPPER :Dict[str,str] = {
    '1': 'HOH',
    '4': 'ATP',
    '22': 'PEO',
    '24': 'SAM',
    '36': 'SAH',
    '1029': 'NIO'
}
"""TODO(CJ)"""

# _DF = pickle.load(open('/panfs/accrepfs.vampire/home/jurichc/EnzyHTP/enzy_htp/chemical/db.pickle', 'rb'))
_DF = {}
"""TODO(CJ)"""

def map_to_pdb(ligand_dict) -> str:

    mask = _DF.InChIKey==ligand_dict['InChIKey']

    if sum(mask):
        assert sum(mask) == 1
        return _DF[mask].iloc[0]['code']

    mask = _DF['name'] == ligand_dict['name']

    if sum(mask):
        print('here')


    mask = _DF['name'] == ligand_dict['name'].upper()

    if sum(mask):
        if sum(mask) == 1:
            return _DF[mask].iloc[0]['code']
        else:
            print(_DF[mask])

    mask = _DF['name'] == ligand_dict['BRENDA Name'].upper()

    if sum(mask):
        if sum(mask) == 1:
            return _DF[mask].iloc[0]['code']
        else:
            print(_DF[mask])


    mask = _DF['name'] == ligand_dict['BRENDA Name']

    if sum(mask):
        if sum(mask) == 1:
            return _DF[mask].iloc[0]['code']
        else:
            print(_DF[mask])



    mask = _DF['formula'] == ligand_dict['formula']

    if sum(mask):
        if sum(mask) == 1:
            return _DF[mask].iloc[0]['code']
        else:
            sliced = _DF[mask].reset_index(drop=True)
            max_tol=len(ligand_dict['BRENDA Name'])*0.33
            distances = list()
            for i, row in sliced.iterrows():
                distances.append(editdistance.eval(row['name'].upper(), ligand_dict['BRENDA Name'].upper()))

            #TODO(CJ): add some kind of filter
            if np.min(distances) <= max_tol:
                return sliced.iloc[np.argmin(distances)]['code']

    return 'unmatched'


def valid_ec_code( code : str ) -> bool:
    """Checks if the supplied enzyme commission (EC) code is valid. Check is simplistic and
    does not consider the supplied code is real. Must be four integers separated by three dots."""


    if code.count('.') != 3:
        return False

    for token in code.split('.'):

        if not token.isnumeric():
            #TODO(CJ): log the error here
            return False

        if not int(token):
            #TODO(CJ): log the error here
            return False

    return True


def convert_formula( raw : str ) -> str:

    temp = [raw[0]]

    for last, curr in zip(raw[:-1], raw[1:]):
        if (curr.isalpha() and last.isnumeric()) or (
            (curr.isalpha() and curr.isupper()) and (last.isalpha() and last.isupper())
            ):
            temp.append(' ')

        temp.append(curr)

    return ''.join(temp)




def parse_ligand(tag) -> Dict:
    """
    """
    result = {}
    child = list(tag.children)[0]
    id_number=child['href'].split('=')[-1]
    
    ligand_url:str=f"https://www.brenda-enzymes.org/molfile.php?LigandID={id_number}"
    r = requests.get(ligand_url)
    if not r.ok:
        #TODO(CJ): error code here if it doesn't load
        pass
    
    result['smiles'] = AllChem.MolToSmiles(AllChem.MolFromMolBlock(r.text, sanitize=False, removeHs=False))


    result['name'] = child.string
    id_map = BRENDA_MAPPER.get(id_number, None)
    if id_map is not None:
        _LOGGER.info(f"Able to map BRENDA id {id_number} to PDB residue code {id_map}!")
        result['pdb'] = id_map
        return result

    url:str=f"https://www.brenda-enzymes.org/ligand.php?brenda_ligand_id={id_number}"
    structure_download:str=f"https://www.brenda-enzymes.org/molfile.php?LigandID={id_number}"

    r = requests.get(url)

    if not r.ok:
        #TODO(CJ): some kind of better error goes here
        assert False

    bs = BeautifulSoup(r.text, features='lxml')

    tokens = list()
    for bb in bs.find_all('div', 'equal'):
        for row in bb.children:
            for child in row.children:
                tokens.append(child.text)
        break


    len_tokens = len(tokens)
    assert len_tokens%2 == 0

    mapper = dict(zip(tokens[:int(len_tokens/2)],tokens[int(len_tokens/2):]))
    result.update(mapper)

    result['formula'] = convert_formula(result['Molecular Formula'])

    result['pdb'] = map_to_pdb( result )

    #TODO(CJ): here is where we map to the ligand

    return result


def parse_ec_number( code : str ) -> Dict:
    """Given a formatted enzyme commission (E.C.) number, 

    Args:

    Returns:
        
    """

    if not valid_ec_code( code ):
        pass
        #TODO(CJ): some error stuff here


    url:str=f"https://www.brenda-enzymes.org/enzyme.php?ecno={code}"


    r = requests.get(url)

    if not r.ok:
        _LOGGER.error(f"The supplied enzyme commission code {code} is invalid in its current format. Exiting...")
        exit( 1 )

    bs = BeautifulSoup(r.text, features='lxml')

    for tag in bs.find_all('div', 'reactcontainer'):
        is_reactant = True
        reactants, products = [], []
        for re_eq in tag.find_all('div', 'react_equation'):
            for child in re_eq.children:
                if child.string and (not child.string.strip()):
                    continue

                if child['class'][0] == 'compound_name':
                    if is_reactant:
                        reactants.append(parse_ligand(child))
                    else:
                        products.append(parse_ligand(child))

                if child['class'][0] == 'react_sign':
                    if child.string == '=':
                        is_reactant = False

    reaction = {
        'products':products,
        'reactants':reactants,
        'ec':code
    }

    return reaction
