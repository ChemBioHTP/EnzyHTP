# Plan 3: Alternative PDB I/O Strategy (Minimal Impact)

## Problem Analysis  
The deepcopy issue only occurs in PDB I/O operations during the `get_file_str()` method. Rather than fixing the deepcopy issue itself, we can modify the PDB I/O to avoid needing deepcopy entirely.

## Solution Overview
Rewrite PDB I/O logic to work directly with the original structure without requiring deepcopy, eliminating the root cause of the problem.

## Implementation Steps

### 1. Eliminate Deepcopy in PDB I/O
- **File**: `enzy_htp/structure/structure_io/pdb_io.py`
- **Action**: Modify `get_file_str()` method:
  - Remove `stru = copy.deepcopy(stru)` line
  - Work directly with original structure  
  - Use temporary modifications with rollback instead of copying

### 2. Implement Temporary Modification Pattern
- **File**: `enzy_htp/structure/structure_io/pdb_io.py`
- **Action**: Add utility functions:
  - `with_temporary_renumbering(structure, func)` - context manager
  - `with_temporary_atom_names(structure, func)` - context manager
  - Apply modifications, run operation, restore original state

### 3. Add Structure State Backup/Restore
- **File**: `enzy_htp/structure/structure.py`
- **Action**: Add methods:
  - `backup_state()` -> returns state dict with atom numbers, names, etc.
  - `restore_state(state_dict)` -> restores structure to previous state
  - `create_state_snapshot()` -> lightweight state capture

### 4. Update Related I/O Operations
- **File**: `enzy_htp/structure/structure_io/pdb_io.py`
- **Action**: Review and update other methods that use deepcopy:
  - `save_structure()` - use temporary modification pattern
  - Any other I/O methods that copy structures
  - Ensure consistent approach across all I/O

### 5. Add Safety Checks
- **File**: `enzy_htp/structure/structure_io/pdb_io.py`
- **Action**: Add validation:
  - Verify structure state is restored after operations
  - Add warnings if structure modification is detected
  - Include rollback mechanisms for failure cases

## Pros  
- ✅ Minimal code changes required
- ✅ Fixes the immediate problem quickly
- ✅ No impact on connectivity system
- ✅ Preserves all existing functionality
- ✅ Low risk of introducing bugs

## Cons
- ❌ Doesn't solve the underlying deepcopy issue
- ❌ May make PDB I/O code more complex
- ❌ Could impact performance with state backup/restore
- ❌ Other parts of codebase may still hit deepcopy issues

## Risk Assessment: Low
- Changes are localized to PDB I/O
- Easy to test and validate
- Can be implemented quickly
- Fallback option if other solutions fail

## Timeline: 1 week