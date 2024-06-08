from pathlib import Path
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem as _rchem
from rdkit.Chem import rdMolTransforms as rdmt

from enzy_htp import config, interface


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

def mimic_torsions(template, reactant):
    template_path = template
    ligand = AllChem.MolFromMol2File(reactant, removeHs=False,  cleanupSubstructures=False)
    template = AllChem.MolFromMol2File(template, removeHs=False,  cleanupSubstructures=False)


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
            assert False


    for (a1,a2,a3,a4) in torsions:
        m1,m2,m3,m4 = template_to_ligand[a1],template_to_ligand[a2],template_to_ligand[a3],template_to_ligand[a4]
        rdmt.SetDihedralDeg(ligand.GetConformer(), m1, m2, m3, m4,
            rdmt.GetDihedralDeg(template.GetConformer(), a1, a2, a3, a4)
        )




    session = interface.pymol.new_session()
    args = [('load', reactant)]
    for aidx, atom in enumerate(ligand.GetAtoms()):
        pos = ligand.GetConformer().GetAtomPosition(aidx)
        args.append(('alter_state', 1, f'rank {aidx}', f"(x,y,z) = ({pos.x},{pos.y},{pos.z})"))

    args.extend([
    ('load', template_path),
    ('align', Path(reactant).stem, Path(template_path).stem),
    ('delete', Path(template_path).stem),
    ('save', reactant)
    ])    

    interface.pymol.general_cmd(session, args)

    return reactant
