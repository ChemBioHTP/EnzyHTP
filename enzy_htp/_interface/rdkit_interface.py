"""
Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2023-10-09
"""

#TODO(CJ): maybe wrap this in a try,except loop
import importlib
from copy import deepcopy
#import rdkit
#from rdkit import Chem as _rchem
#import rdkit.Chem.AllChem as _rchem_ac
from typing import List
from pathlib import Path

from .base_interface import BaseInterface

from enzy_htp import _LOGGER

from enzy_htp.core.general import HiddenPrints
from enzy_htp.core import file_system as fs
from enzy_htp._config.rdkit_config import RDKitConfig, default_rdkit_config

from enzy_htp.structure import (
    Ligand,
    Mol2Parser
)

class RDKitInterface(BaseInterface):
    pass

    def __init__(self, parent, config: RDKitConfig = None) -> None:
        """Simplicstic constructor that requires the parent interface as an argument and optionally takes an RDKitConfig instance.
        Calls parent constructor.
        """
        super().__init__(parent, config, default_rdkit_config)
        self.chem_ = None
        try:
            self.chem_ = importlib.import_module('rdkit.Chem')
        except:
            pass

        from rdkit import RDLogger #TODO(CJ): update this
        RDLogger.DisableLog('rdApp.*')

    @property
    def chem(self) -> "module":
        """Gets the modeller module if it exists, raises an error if not."""
        if self.chem_ is None:
            err_msg:str="The 'rdkit.Chem' python package is not installed in this environment. Cannot use RDKitInterface()."
            _LOGGER.error( err_msg )  
            raise ImportError(err_msg)
        return self.chem_

    def _supported_ftype(self, molfile: str) -> None:
        """Is the supplied file type supported for use in rdkit? Errors and exists if not."""
        ext: str = Path(molfile).suffix

        if ext not in self.config_.SUPPORTED_FTYPES:
            _LOGGER.error(
                f"The supplied file {molfile} has an unsupported extention '{ext}'. RDKitInterface supports: {', '.join(self.config_.SUPPORTED_FTYPES)}. Exiting..."
            )
            exit(1)

    def check_rdkit_installed(self) -> None:
        """ """
        if "rdkit" in self.missing_py_modules():
            _LOGGER.error("rdkit is NOT installed. Use 'conda install -c conda-forge -y -q rdkit'")
            exit(1)

    def _load_molecule(self,
                       molfile: str,
                       sanitize: bool = False,
                       removeHs: bool = False,
                       cleanupSubstructures: bool = False,
                       proximityBonding: bool = False,
                       strictParsing: bool = False) :
        """ """
        #TODO(CJ): add in checks that the file exists
        #TODO(CJ): get the list of molecules out if needed
        fs.check_file_exists(molfile)
        self._supported_ftype(molfile)

        ext: str = Path(molfile).suffix

        if ext == '.mol2':
            return _rchem.MolFromMol2File(molfile, sanitize, removeHs, cleanupSubstructures)

        if ext == '.pdb':
            return _rchem.Mol2FromPDBFile(molfile, sanitize, removeHs, 0, proximity_bonding)

        if ext == '.mol':
            return _rchem.MolFromMolFile(molfile, sanitize, removeHs, strictParsing)

        if ext == '.sdf':
            reader = _rchem.SDMolSupplier(molfile)
            for rr in reader:
                return rr

    def _save_molecule(self, mol, outfile: str, kekulize: bool = True) -> str:
        """ """
        self._supported_ftype(outfile)

        ext: str = Path(outfile).suffix

        if ext == '.mol' or ext == '.mol2':
            temp = str(Path(outfile).with_suffix('.mol'))
            _rchem.MolToMolFile(mol, temp, kekulize=kekulize)
            if ext == '.mol2':
                session = self.parent.pymol.new_session()
                self.parent.pymol.convert(session, temp, new_ext='.mol2')

        if ext == '.pdb':
            _rchem.MolToPDBFile(mol, outfile)

        if ext == '.sdf':
            writer = _rchem.SDWriter(outfile)
            writer.SetKekulize(kekulize)
            writer.write(mol)
            writer.close()

        return outfile

    def mol_from_ligand(self, 
                ligand : Ligand, 
                removeHs:bool=True,
                cleanupSubstructures:bool=True,
                work_dir:str=None ) -> "rdkit.Chem.Mol":
        if not work_dir:
            work_dir = self.parent.config['system.SCRATCH_DIR']
   
        fs.safe_mkdir( work_dir )

        temp_file:str = f"{work_dir}/rdkit_temp.mol2"
        
        lp = Mol2Parser()

        lp.save_ligand( temp_file, ligand )

        result = self.chem.MolFromMol2File( temp_file, removeHs=removeHs, cleanupSubstructures=cleanupSubstructures )

        fs.safe_rm( temp_file )

        return result

    def find_mcs(self, ligand1:Ligand, ligand2:Ligand):


        mol1 = self.mol_from_ligand( ligand1 )
        mol2 = self.mol_from_ligand( ligand2 )
       
        module = None
        try:
            module = importlib.import_module('rdkit.Chem.rdFMCS')
        except:
            pass

        assert module

        return  module.FindMCS( [mol1, mol2] )


    def kekulize(self, molfile: str, outfile: str) -> str:
        """ """
        mol = self._load_molecule(molfile)
        #TODO(CJ): need to add some stuff for phosphates/carboxylic acids
        self._save_molecule(mol, outfile)
        return outfile

    def num_rotatable_bonds(self, molfile: str) -> int:
        """How many rotatble bonds does the molecule have?"""
        mol = self._load_molecule(molfile)

        return _rchem.rdMolDescriptors.CalcNumRotatableBonds(mol)

    def volume(self, molfile: str, gridSpacing: float = 0.2, boxMargin: float = 2.0) -> float:
        """Calculates the volume of the molecule in the supplied file. All calculations done in Angstroms.

        Args:
            molfile:
            gridSpacing:
            boxMargin:

        Returns:
            The volume in Angstroms^3.
            
        """
        self.check_rdkit_installed()
        mol: _rdkit.Chem = self._load_molecule(molfile)
        return _rchem_ac.ComputeMolVolume(mol, gridSpacing=gridSpacing, boxMargin=boxMargin)

    def update_ligand_positions(self, ligand:Ligand, mol, cidx:int=-1) -> None:
        
        assert len(ligand.atoms) == len(mol.GetAtoms())
       
        for idx, latom in enumerate(ligand.atoms):
            pos = mol.GetConformer(cidx).GetAtomPosition(idx)
            latom.coord = (pos.x, pos.y, pos.z)


    def generate_conformers( self,
            ligand:Ligand, 
            n_conformers:int,
            rms_cutoff:float,
            attempts:int,
            rng:int=1996 ) -> List[Ligand]:

        module = None
        try:
            module = importlib.import_module('rdkit.Chem.rdDistGeom')
        except:
            pass

        assert module
        mol = self.mol_from_ligand(ligand, removeHs=False )
        
        module.EmbedMultipleConfs( 
                mol,
                numConfs=n_conformers,
                maxAttempts=attempts,
                randomSeed=rng,
                pruneRmsThresh=rms_cutoff,
                )
        
        result = list()
        for idx in range( mol.GetNumConformers() ):
            temp = deepcopy( ligand )
            self.update_ligand_positions( temp, mol, idx )
            result.append( temp )

        return result

    def apply_coords(self,
                template:Ligand,
                ligand:Ligand) -> None:
        from rdkit.Chem import AllChem
        template_to_ligand = dict()
        atom_names = list()
        tmol = self.mol_from_ligand(template, False)
        lmol = self.mol_from_ligand(ligand, False)

        orig_smiles=AllChem.MolToSmiles(lmol,canonical=True)
        for aidx, at in enumerate(tmol.GetAtoms()):
            target_name:str=at.GetPropsAsDict()['_TriposAtomName']
            atom_names.append(target_name)
            for lidx, al in enumerate(lmol.GetAtoms()):
                if target_name == al.GetPropsAsDict()['_TriposAtomName']:
                    template_to_ligand[aidx] = lidx
                    break
            else:
                #TODO(CJ): put an error code here
                pass
                #assert False, target_name
        
        lconf = lmol.GetConformer()
        tconf = tmol.GetConformer()
        
        ff = AllChem.UFFGetMoleculeForceField(lmol)
        for tidx, lidx in template_to_ligand.items():
            
            lconf.SetAtomPosition(lidx, tconf.GetAtomPosition(tidx))

            ff.UFFAddPositionConstraint(lidx, 0.05, 10000)
        ff.Minimize()
        from rdkit.Chem import rdMolTransforms as rdmt

        #ff = AllChem.UFFGetMoleculeForceField(lmol)
        for idx in range(lmol.GetNumAtoms()):
            atom = lmol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 1:
                ff.UFFAddPositionConstraint(idx, 0.05, 100)
            else:
                pidx = None
                for bb in lmol.GetBonds():
                    bidxs = [bb.GetBeginAtomIdx(), bb.GetEndAtomIdx()]
                    if idx == bidxs[0]:
                        pidx = bidxs[1]
                    
                    if idx == bidxs[1]:
                        pidx = bidxs[0]
                

                rdmt.SetBondLength(lconf, pidx, idx, 1.10)
               
        ff.Minimize()
        self.update_ligand_positions(ligand, lmol)
        end_smiles=AllChem.MolToSmiles(lmol, canonical=True)
        assert orig_smiles == end_smiles, f"{orig_smiles}\n{end_smiles}"
        #assert False
