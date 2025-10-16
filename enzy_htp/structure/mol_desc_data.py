"""A dataclass to hold generic molecular description data for file format conversion."""
from dataclasses import dataclass, field, fields
from typing import Dict, List, Any


@dataclass
class MolDescData:
    """A generic data structure for molecular descriptions.

    This class serves as an intermediate representation for converting between
    different molecular file formats. It stores essential information about
    atoms, bonds, and the overall molecule.

    Attributes:
        name (str): The name of the molecule.
        atoms (List[Dict]): A list of dictionaries, each representing an atom.
            Each dictionary should have the following keys:
            - 'id' (int): A unique identifier for the atom.
            - 'atom_name' (str): The name of the atom (e.g., 'CA', 'H1').
            - 'atom_type' (str): The type of the atom (e.g., 'C.3', 'H').
            - 'charge' (float): The partial charge of the atom.
            - 'coords' (List[float]): A list of three floats representing the
              x, y, and z coordinates of the atom.
        bonds (List[Dict]): A list of dictionaries, each representing a bond.
            Each dictionary should have the following keys:
            - 'atom1_id' (int): The ID of the first atom in the bond.
            - 'atom2_id' (int): The ID of the second atom in the bond.
            - 'bond_type' (str): The type of the bond (e.g., 'single', 'aromatic').
    """
    name: str = "MOL"
    charge_type: str = "USER_CHARGES"
    atoms: List[Dict] = field(default_factory=list)
    bonds: List[Dict] = field(default_factory=list)

    def update(self, updates: Dict[str, Any]) -> None:
        """Update the molecular description data in place.

        The method is designed to make it easy for conversion routines to tweak
        metadata (for example, ``charge_type``) or replace the atom/bond collections
        prior to serialisation.

        Args:
            updates: Mapping of field names to their desired values.

        Raises:
            KeyError: If ``updates`` contains a key that is not a valid
                MolDescData attribute.
            TypeError: If a provided value has an unexpected type.
        """
        if not updates:
            return

        valid_fields = {f.name for f in fields(self)}

        for key, value in updates.items():
            if key not in valid_fields:
                raise KeyError(f"MolDescData does not support field '{key}'.")

            if key in {"atoms", "bonds"}:
                if value is None:
                    continue
                if not isinstance(value, list):
                    raise TypeError(
                        f"MolDescData field '{key}' expects a list, received {type(value).__name__}."
                    )
                setattr(self, key, list(value))
                continue

            setattr(self, key, value)
