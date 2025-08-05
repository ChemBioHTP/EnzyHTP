# Plan 4: Hybrid Approach with Connectivity Modes

## Problem Analysis
The connectivity system is essential for many operations but causes deepcopy issues. We need a solution that provides connectivity when needed but allows safe copying when required.

## Solution Overview
Implement a dual-mode connectivity system that can operate in both "connected" and "copyable" modes, with seamless transitions between them.

## Implementation Steps

### 1. Add Connectivity Mode Enum
- **File**: `enzy_htp/structure/structure_enchantment/connectivity.py`
- **Action**: Define connectivity modes:
```python
class ConnectivityMode(Enum):
    DISCONNECTED = "disconnected"  # No connections, safe for deepcopy
    CONNECTED = "connected"        # Full connections, not safe for deepcopy  
    SERIALIZED = "serialized"      # Connections stored as data, safe for deepcopy
```

### 2. Implement Mode Switching
- **File**: `enzy_htp/structure/structure.py`
- **Action**: Add methods:
  - `set_connectivity_mode(mode: ConnectivityMode)`
  - `get_connectivity_mode() -> ConnectivityMode`  
  - `ensure_mode(required_mode: ConnectivityMode)`
  - Auto-convert between modes as needed

### 3. Serialized Connection Storage
- **File**: `enzy_htp/structure/atom.py`
- **Action**: Add connection serialization:
  - `_connection_data: Dict` - stores connection info as data
  - `_live_connections: List` - actual circular references
  - `serialize_connections()` - convert refs to data
  - `deserialize_connections()` - convert data to refs

### 4. Context Managers for Safe Operations
- **File**: `enzy_htp/structure/structure_enchantment/connectivity.py`
- **Action**: Add context managers:
```python
@contextmanager
def copyable_mode(structure):
    """Temporarily switch to copyable mode"""
    original_mode = structure.get_connectivity_mode()
    structure.set_connectivity_mode(ConnectivityMode.SERIALIZED)
    try:
        yield structure
    finally:
        structure.set_connectivity_mode(original_mode)
```

### 5. Update PDB I/O with Context Manager
- **File**: `enzy_htp/structure/structure_io/pdb_io.py`
- **Action**: Modify deepcopy operations:
```python
def get_file_str(stru, if_renumber, if_fix_atomname):
    with copyable_mode(stru):
        stru_copy = copy.deepcopy(stru)
        # ... rest of operation
```

### 6. Update Connectivity-Dependent Code
- **File**: Various files using connectivity
- **Action**: Update code to ensure correct mode:
  - Parameterization code - ensure connected mode
  - Analysis code - ensure connected mode  
  - I/O operations - use copyable mode context

## Pros
- ✅ Best of both worlds - connectivity when needed, copying when needed
- ✅ Backwards compatible - existing code works unchanged
- ✅ Flexible - can optimize for different use cases
- ✅ Clear mode transitions with context managers
- ✅ Fixes deepcopy issue comprehensively

## Cons
- ❌ Complex implementation with multiple modes
- ❌ Potential for mode confusion and bugs
- ❌ Performance overhead from mode switching
- ❌ Requires updates to connectivity-dependent code
- ❌ More complex mental model for developers

## Risk Assessment: Medium-High
- Complex system with multiple states
- Could introduce subtle bugs if mode switching is wrong
- Requires thorough testing of all mode combinations
- Performance implications need careful consideration

## Timeline: 4-5 weeks