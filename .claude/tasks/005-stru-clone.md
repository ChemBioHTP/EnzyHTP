# Structure Clone Method Implementation Plan

## Problem Analysis

The current clone methods in Structure, Chain, and Residue classes create fast copies without using deepcopy, but the connectivity cloning is incomplete. The TODO comments indicate that connectivity should be handled by:

1. Creating a mapping between old atoms and new atoms
2. For each atom in the new structure, cloning connectivity from the corresponding old atom
3. Warning and ignoring connections to atoms not present in the new structure

## Current State

- **Atom.clone()**: Already complete - creates atoms without connectivity
- **Residue.clone()**: Partially implemented - handles intra-residue connectivity only
- **Chain.clone()**: Has TODO placeholder for connectivity handling
- **Structure.clone()**: Has TODO placeholder for connectivity handling

## Connectivity Structure

Based on the code analysis:
- `Atom.connect`: List of tuples `[(connected_atom, bond_type), ...]`
- Bond types are strings (e.g., "single", "double") or None
- Connectivity is bidirectional (both atoms reference each other)

## Implementation Plan

### Phase 1: Core Mapping Infrastructure
1. **Create atom mapping methods** at each level:
   - `Structure().create_atom_mapping(self, other: Structure) -> Dict[Atom, Atom]`
   - `Chain().create_atom_mapping(self, other: Chain) -> Dict[Atom, Atom]` 
   - `Residue().create_atom_mapping(self, other: Residue) -> Dict[Atom, Atom]`

### Phase 2: Connectivity Cloning Implementation
2. **Fix Residue.clone()**: 
   - Remove existing intra-residue connectivity code (it's buggy)
   - Use mapping-based approach for consistency

3. **Implement Chain.clone()**: 
   - Create atom mapping between old and new chains
   - Clone connectivity for atoms within the chain scope
   - Handle inter-residue bonds within chain

4. **Implement Structure.clone()**:
   - Create atom mapping between old and new structures  
   - Clone connectivity for all atoms in structure scope
   - Handle inter-chain bonds
   - Warn about missing connections to atoms outside structure

### Phase 3: Testing & Validation
5. **Write comprehensive tests**:
   - Test intra-residue connectivity preservation and the peptide bond connection to other residues is not.
   - Test inter-residue connectivity within chains
   - Test partial structure cloning (subset of original)
   - Test warning behavior for out-of-scope connections

## Technical Details

### Atom Mapping Strategy
- Make `is_topology_subset_atomic` method and use it to do a san check before creating the map. (this method mimic `is_same_topology_atomic` but use set and issubset) 
- Use atom keys (`chain.residue_idx.atom_name`) for mapping
- Handle cases where new structure is a subset of old structure
- Provide clear error messages for mapping failures

### Connectivity Cloning Algorithm
The connectivity cloning will be implemented as a class method in the Atom class:

```python
@classmethod 
def clone_connectivity(cls, atom_mapping: Dict[Atom, Atom]) -> None:
    """Clone connectivity between atoms using the provided mapping.
    
    Args:
        atom_mapping: Dictionary mapping old atoms to new atoms
    """
    for old_atom, new_atom in atom_mapping.items():
        if not old_atom.is_connected():
            continue
            
        new_connections = []
        for connected_old_atom, bond_type in old_atom.connect:
            if connected_old_atom in atom_mapping:
                connected_new_atom = atom_mapping[connected_old_atom] 
                new_connections.append((connected_new_atom, bond_type))
            else:
                # Log debug info about missing connection (expected for partial cloning)
                _LOGGER.debug(f"Skipping connection from {old_atom.key} to {connected_old_atom.key} - target not in cloned structure")
        
        new_atom.connect = new_connections if new_connections else None
```

### Error Handling
- Warn when connections reference atoms not in the new structure
- Provide informative error messages with atom keys
- Ensure bidirectional connectivity is maintained

## Implementation Order
1. Add `is_topology_subset_atomic()` method to Structure class
2. Implement `Atom.clone_connectivity()` class method
3. Implement `create_atom_mapping()` methods:
   - Structure.create_atom_mapping(other: Structure)
   - Chain.create_atom_mapping(other: Chain) 
   - Residue.create_atom_mapping(other: Residue)
4. Fix Residue.clone() connectivity (remove buggy code, use mapping + Atom.clone_connectivity())
5. Implement Chain.clone() connectivity handling
6. Implement Structure.clone() connectivity handling
7. Write and run comprehensive tests

This approach ensures consistency across all levels and handles the complex case of partial structure cloning while maintaining clear separation of concerns.