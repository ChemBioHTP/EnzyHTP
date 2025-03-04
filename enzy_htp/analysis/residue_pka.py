from typing import List
from enzy_htp.core.logger import _LOGGER

# Importing PROPKA components
from propka.lib import loadOptions
from propka.input import read_parameter_file, read_molecule_file
from propka.parameters import Parameters
from propka.molecular_container import MolecularContainer


def get_residue_pka(target_mask: str, pdb_path: str) -> List[float]:
    try:
        # Parse target residues
        target_resis = list(map(int, target_mask.removeprefix(':').split(',')))

        # Load PROPKA options
        options = loadOptions([pdb_path])  # Defaults for simplicity
        parameters = read_parameter_file(options.parameters, Parameters())
        my_molecule = MolecularContainer(parameters, options)

        # Read and initialize molecule
        my_molecule = read_molecule_file(pdb_path, my_molecule)

        # Calculate pKa values
        my_molecule.calculate_pka()
        conformation = my_molecule.conformations.get('AVR', {})
        residues_pka = {}

        # Extract pKa data for all residues
        for group in conformation.groups:
            atom = group.atom
            residues_pka[atom.res_num] = group.pka_value

        # Return pKa values for target residues
        return [residues_pka[res] for res in target_resis if res in residues_pka]

    except FileNotFoundError:
        _LOGGER.error(f"PDB file not found: {pdb_path}")
        raise
    except Exception as e:
        _LOGGER.error(f"Error in get_residue_pka: {e}")
        raise

