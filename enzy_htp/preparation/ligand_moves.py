"""Implements some basic ligand manipulation and movement functionality.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2024-06-10
"""
from pathlib import Path
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem as _rchem
from rdkit.Chem import rdMolTransforms as rdmt

from rdkit.Chem.rdMolAlign import AlignMol

from enzy_htp import config, interface

from enzy_htp.structure import Ligand


def enumerate_torsions(mol):
    ri = mol.GetRingInfo()

    ring_mapper = defaultdict(list)
    for ridx, ring in enumerate(ri.AtomRings()):
        ring_name = f"ring_{ridx}"
        for rr in ring:
            ring_mapper[rr].append(ring_name)
    torsionSmarts = '[!$(*#*)&!D1]~[!$(*#*)&!D1]'
    torsionQuery = Chem.MolFromSmarts(torsionSmarts)
    matches = mol.GetSubstructMatches(torsionQuery)
    torsionList = []
    for match in matches:
        idx2 = match[0]
        idx3 = match[1]
        bond = mol.GetBondBetweenAtoms(idx2, idx3)
        jAtom = mol.GetAtomWithIdx(idx2)
        kAtom = mol.GetAtomWithIdx(idx3)
        if (((jAtom.GetHybridization() != Chem.HybridizationType.SP2)
            and (jAtom.GetHybridization() != Chem.HybridizationType.SP3))
            or ((kAtom.GetHybridization() != Chem.HybridizationType.SP2)
            and (kAtom.GetHybridization() != Chem.HybridizationType.SP3))):
            continue
        for b1 in jAtom.GetBonds():
            if (b1.GetIdx() == bond.GetIdx()):
                continue
            idx1 = b1.GetOtherAtomIdx(idx2)
            for b2 in kAtom.GetBonds():
                if ((b2.GetIdx() == bond.GetIdx())
                    or (b2.GetIdx() == b1.GetIdx())):
                    continue
                idx4 = b2.GetOtherAtomIdx(idx3)
                # skip 3-membered rings
                if (idx4 == idx1):
                    continue
                if idx2 in ring_mapper and idx3 in ring_mapper  and set(ring_mapper[idx2])&set(ring_mapper[idx3]):
                    continue
                torsionList.append((idx1, idx2, idx3, idx4))
    return torsionList

def mimic_torsions(t_ligand:Ligand, r_ligand:Ligand) -> None:
    
    template = interface.rdkit.mol_from_ligand(t_ligand, removeHs=False, cleanupSubstructures=False)
    ligand  = interface.rdkit.mol_from_ligand(r_ligand, removeHs=False, cleanupSubstructures=False)

    torsions = enumerate_torsions(template)

    template_to_ligand = dict()

    atom_names = list()
    for aidx, at in enumerate(template.GetAtoms()):
        target_name:str=at.GetPropsAsDict()['_TriposAtomName']
        atom_names.append(target_name)
        for lidx, al in enumerate(ligand.GetAtoms()):
            if target_name == al.GetPropsAsDict()['_TriposAtomName']:
                template_to_ligand[aidx] = lidx
                break
        else:
            #TODO(CJ): put an error code here
            pass
            #assert False, target_name


    for (a1,a2,a3,a4) in torsions:
        if not (a1 in template_to_ligand and a2 in template_to_ligand and a3 in template_to_ligand and a4 in template_to_ligand):
            continue
        m1,m2,m3,m4 = template_to_ligand[a1],template_to_ligand[a2],template_to_ligand[a3],template_to_ligand[a4]
        rdmt.SetDihedralDeg(ligand.GetConformer(), m1, m2, m3, m4,
            rdmt.GetDihedralDeg(template.GetConformer(), a1, a2, a3, a4)
        )

    amapper = list(zip(template_to_ligand.values(), template_to_ligand.keys()))
    AlignMol( ligand, template, atomMap=amapper ) #TODO(CJ): put into the rdkit interface

    interface.rdkit.update_ligand_positions(r_ligand, ligand)


def ligand_mcs_score( l1, l2 ) -> float:
    
    if l1 is None or l2 is None:
        return 0.0

    ct = 0
    for aa in l1.atoms:
        if aa.element == 'H':
            continue
        ct += 1
    
    for aa in l2.atoms:
        if aa.element == 'H':
            continue
        ct += 1
    mcs_result = interface.rdkit.find_mcs( l1 , l2 )
    return mcs_result.numAtoms /  (ct - mcs_result.numAtoms)


def mimic_torsions_mcs(t_ligand:Ligand, r_ligand:Ligand) -> None:
    
    mcs = interface.rdkit.find_mcs( t_ligand, r_ligand )

    template = interface.rdkit.mol_from_ligand(t_ligand, removeHs=False, cleanupSubstructures=False)
    ligand  = interface.rdkit.mol_from_ligand(r_ligand, removeHs=False, cleanupSubstructures=False)

    ss1 = template.GetSubstructMatch(mcs.queryMol)
    ss2 = ligand.GetSubstructMatch(mcs.queryMol)


    template_to_ligand = dict()
    
    template_to_ligand = dict( zip (
            ss1, ss2
    ))

    torsions = enumerate_torsions(template)
    for (a1,a2,a3,a4) in torsions:
        if a1 not in template_to_ligand:
            continue
        if a2 not in template_to_ligand:
            continue
        if a3 not in template_to_ligand:
            continue
        if a4 not in template_to_ligand:
            continue
        m1,m2,m3,m4 = template_to_ligand[a1],template_to_ligand[a2],template_to_ligand[a3],template_to_ligand[a4]
        rdmt.SetDihedralDeg(ligand.GetConformer(), m1, m2, m3, m4,
            rdmt.GetDihedralDeg(template.GetConformer(), a1, a2, a3, a4)
        )

    amapper = list(zip(template_to_ligand.values(), template_to_ligand.keys()))
    AlignMol( ligand, template, atomMap=amapper ) #TODO(CJ): put into the rdkit interface
    

    interface.rdkit.update_ligand_positions(r_ligand, ligand)

