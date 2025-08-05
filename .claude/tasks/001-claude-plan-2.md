# Plan 2: Lazy Connectivity Initialization

## Problem Analysis
The connectivity system creates circular references immediately upon calling `init_connectivity()`, which then prevents deepcopy operations. The issue is that connectivity is initialized too early in the workflow.

## Solution Overview
Modify the connectivity system to use lazy initialization, only creating circular references when connectivity data is actually needed.

## Implementation Steps

### 1. Add Connectivity State Management
- **File**: `enzy_htp/structure/structure.py`
- **Action**: Add structure-level connectivity state:
  - `_connectivity_initialized: bool = False`
  - `_connectivity_data: Dict = None` (cached connectivity info)
  - `connectivity_mode: str = "lazy"` (lazy vs eager)

### 2. Modify Connectivity Initialization
- **File**: `enzy_htp/structure/structure_enchantment/connectivity.py`
- **Action**: Change `init_connectivity()` to:
  - Store connectivity calculation results without creating circular refs
  - Set `_connectivity_initialized = True`
  - Add `_ensure_connectivity()` method for lazy loading

### 3. Lazy Connection Properties
- **File**: `enzy_htp/structure/atom.py`
- **Action**: Modify connection-related methods:
  - `is_connected()` -> check cached data, create refs if needed
  - `get_connections()` -> lazy load and return connections
  - `connect_to()` -> ensure connectivity before connecting

### 4. Update PDB I/O Workflow
- **File**: `enzy_htp/_interface/amber_interface.py`
- **Action**: Modify `_write_combining_tleap_input`:
  - Call `structure.prepare_for_copy()` before deepcopy
  - Use `structure.restore_connectivity()` after deepcopy
  - Add temporary connectivity suspension

### 5. Add Copy-Safe Mode
- **File**: `enzy_htp/structure/structure.py` 
- **Action**: Add methods:
  - `suspend_connectivity()` - temporarily breaks circular refs
  - `restore_connectivity()` - rebuilds connections from cached data
  - `prepare_for_copy()` - convenience method for safe copying

## Pros
- ✅ Minimal changes to existing code
- ✅ Preserves performance when connectivity not needed
- ✅ Clean separation of concerns
- ✅ Fixes issue without complex deepcopy logic

## Cons
- ❌ Changes fundamental behavior of connectivity system
- ❌ May cause subtle bugs if connectivity checks are wrong
- ❌ Could impact code that assumes immediate connectivity
- ❌ Requires updates to all connectivity-dependent code

## Risk Assessment: High
- Changes core behavior that other code may depend on
- Could introduce hard-to-debug timing issues
- Requires extensive testing of connectivity-dependent features

## Timeline: 3-4 weeks