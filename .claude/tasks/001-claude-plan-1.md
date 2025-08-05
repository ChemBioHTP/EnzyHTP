# Plan 1: Implement Custom Deepcopy for Structure Classes

## Problem Analysis
The connectivity system in EnzyHTP creates circular references between atoms that cause Python's `copy.deepcopy()` to fail with a RecursionError. This prevents PDB I/O operations that rely on deepcopy for structure manipulation.

## Solution Overview  
Implement custom `__deepcopy__` methods for core structure classes to handle circular references gracefully.

## Implementation Steps

### 1. Add Custom Deepcopy to Base Classes
- **File**: `enzy_htp/structure/atom.py`
- **Action**: Add `__deepcopy__` method to `Atom` class that:
  - Creates new atom with same properties
  - Temporarily breaks circular references during copy
  - Restores connections after all atoms are copied

### 2. Add Custom Deepcopy to Structure Class
- **File**: `enzy_htp/structure/structure.py` 
- **Action**: Add `__deepcopy__` method to `Structure` class that:
  - Creates new structure shell
  - Copies all chains, residues, atoms with custom logic
  - Rebuilds connectivity after copying

### 3. Implement Connectivity-Aware Copying
- **File**: `enzy_htp/structure/structure_enchantment/connectivity.py`
- **Action**: Add utility functions:
  - `break_circular_refs(structure)` - temporarily removes circular refs
  - `restore_connectivity(structure)` - rebuilds connectivity after copy
  - `deepcopy_with_connectivity(structure)` - safe deepcopy wrapper

### 4. Update PDB I/O to Use Custom Deepcopy
- **File**: `enzy_htp/structure/structure_io/pdb_io.py`
- **Action**: Replace `copy.deepcopy()` calls with:
  - `deepcopy_with_connectivity()` for connected structures
  - Regular deepcopy for non-connected structures

## Pros
- ✅ Maintains all existing functionality
- ✅ Fixes deepcopy issues comprehensively  
- ✅ No breaking changes to user code
- ✅ Can be implemented incrementally

## Cons  
- ❌ Requires significant implementation effort
- ❌ Adds complexity to core classes
- ❌ May impact performance of copy operations
- ❌ Needs careful testing to avoid breaking connectivity

## Risk Assessment: Medium
- Implementation complexity is moderate
- Could introduce bugs if connectivity logic is wrong
- Requires deep understanding of circular reference patterns

## Timeline: 2-3 weeks