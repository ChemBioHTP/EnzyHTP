#TODO(CJ) documentation
import pandas as pd
from collections import defaultdict
from biopandas.pdb import PandasPdb
from typing import List, Set, Dict, Tuple


from .atom import Atom
from .metal_atom import MetalAtom
from .residue import Residue
from .solvent import Solvent, residue_to_solvent
from .ligand import Ligand, residue_to_ligand
from .chain import Chain
from .structure import Structure


def structure_from_pdb(fname: str) -> Structure:
    # TODO(CJ) Make this its own file
	# TODO(CJ) check for file validity
    """
    extract the structure from PDB path. Capable with raw experimental and Amber format
    ---------
    input = path (or file or file_str)
        split the file_str to chain and init each chain
    ligand_list: ['NAME',...]
        User specific ligand names. Only extract these if provided. 
    ---------
    Target:
    - structure(w/name)  - chain - residue - atom
                        |- metalatom(atom)
                        |- ligand(residue)
                        |- solvent(residue)
    - ... (add upon usage)
    """
    def name_chains(mapper: defaultdict) -> None:
        def legal_chain_names(mapper) -> List[str]:
            result = list(string.ascii_uppercase)
            taken = set(list(mapper.keys()))
            result = list(filter(lambda s: s not in taken, result))
            return list(reversed(result))

        key_names = set(list(map(lambda kk: kk.strip(), mapper.keys())))
        if "" not in key_names:
            return mapper
        unnamed = list(mapper[""])
        del mapper[""]

        names = legal_chain_names(mapper)
        unnamed = sorted(unnamed, key=lambda r: -r.min_line())
        new_chain: List[Residue] = [unnamed.pop()]

        while len(unnamed):
            new_res = unnamed.pop()
            if new_chain[-1].neighbors(new_res):
                new_chain.append(new_res)
                continue

            new_name = names.pop()
            mapper[new_name] = new_chain
            new_chain = [new_res]

        if len(new_chain):
            new_name = names.pop()
            mapper[new_name] = new_chain

        return mapper

    def build_residues( df: pd.DataFrame ) -> Dict[str,Residue]:
        mapper = defaultdict(list)
        for i, row in df.iterrows():
            aa = Atom(**row)
            mapper[aa.residue_key()].append(aa)
        
        result : Dict[str, Residue] = dict()
        for res_key, atoms in mapper.items():
            result[res_key] = Residue(
                residue_key=res_key, atoms=sorted(atoms, key=lambda a: a.atom_number)
            )
        return result 

    def build_chains(mapper : Dict[str, Residue]) -> Dict[str,Chain]:
        mapper = defaultdict(list)
        for res in mapper.values():
            mapper[res.chain()].append(res)

        mapper = name_chains(mapper)
        result : Dict[str,Chain] = dict()
        # ok this is where we handle missing chain ids
        for chain_name, residues in mapper.items():
            result[chain_name] = Chain(
                chain_name, sorted(residues, key=lambda r: r.num())
            )
        return result
    
    def get_metalatoms(chains : Dict[str,Chain]) -> Tuple[ List[Chain], List[MetalAtom]] :
        metalatoms : List[MetalAtom] = list()
        for cname, chain in chains.items():
            if chain.is_metal():
                metalatoms.append(chain)
                del chains[cname]
        
        if not len(metalatoms):
            return (chains, metalatoms)
        # Break pseudo residues into atoms and convert to Metalatom object
        holders = []
        for pseudo_resi in metalatoms:
            for metal in pseudo_resi:
                holders.append(Metalatom.fromAtom(metal))
        metalatoms = holders

        # clean empty chains
        for i in range(len(raw_chains) - 1, -1, -1):
            if len(raw_chains[i]) == 0:
                del raw_chains[i]

        return raw_chains, metalatoms


    def get_ligands(chains : Dict[str,Chain], ligand_list=None) -> Tuple[Dict[str,Chain],List[Ligand]]: 
        bad_chains = []
        ligands = []
        for cname, chain in chains.items():

            if chain.is_HET():
                continue

            for idx, residue in enumerate(chain.residues()[::-1]):
                if ligand_list is not None and residue.name in ligand_list:
                    _LOGGER.info(
                        f"Structure: Found user assigned ligand in raw {chain.name()} {residue.name} {residue.num_}"
                    )
                    ligands.append(residue)
                    del chain[idx]
                elif not residue.is_rd_solvent() and not residue.is_rd_non_ligand():
                    _LOGGER.warn(
                        f"Structure: Found ligand in raw {chain.name()} {residue.name} {residue.num_}"
                    )
                    ligands.append(residue)
                    del chain[idx]

            if chain.empty():
                bad_chains.append(cname)

        for bc in bad_chains:
            chains[bc]

        return (chains, list(map(residue_to_ligand, ligands)))


    def get_solvents(chains : Dict[str,Chain]) -> Tuple[Dict[str,Chain],List[Solvent]]:
        """
        get solvent from raw chains and clean chains by deleting the solvent part
        -----
        Method: Assume metal/ligand/solvent can anywhere. Base on rd_solvent_list
        """
        solvents = []
        bad_chains = []
        for cname, chain in chains.items():
            for idx, residue in enumerate(chain.residues()[::-1]):
                if residue.is_rd_solvent():
                    _LOGGER.warn(
                        f"Structure: found solvent in raw {residue.name} {residue.id}"
                    )
                    solvents.append(residue)
                    del chain[idx]

            if chain.empty():
                bad_chains.append(cname)
        # Convert pseudo residues to Ligand object
        for bc in bad_chains:
            del chains[bc]

        return ( chains, list(map(residue_to_solvent, solvents)))



    parser = PandasPdb()
    parser.read_pdb(fname)
    res_mapper : Dict[str,Residue] = build_residues( parser.df['ATOM'] )
    chain_mapper : Dict[str, Chain] = build_chains( res_mapper )
    metal_atoms : List[MetalAtom] = list()
    (chain_mapper, metalatoms) = get_metalatoms( chain_mapper )
    (chain_mapper, ligands ) = get_ligands( chain_mapper ) 
    (chain_mapper, solvents ) = get_solvents( chain_mapper )
    result = Structure(
        chains=chain_mapper,
        metal_atoms=chain_mapper,
        solvents=solvents,
        ligands=ligands
	)
    return result
