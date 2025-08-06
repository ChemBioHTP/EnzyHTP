# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Installation and Setup
- Avoid installing the package yourself. If you found a problem involving environment or need to test install, ask for help from the user.
- `./dev-tools/pip-install` - Install enzy_htp from source using pip (do not install dependencies)
- `./dev-tools/conda-install` - Install enzy_htp and dependencies using conda

## Testing

The project uses `pytest` for testing. Use a default timeout wall of 10 min.

**IMPORTANT:** Avoid running `pytest` on the entire test suite, as it is very time-consuming. Instead, run specific tests or smaller groups of tests.

### Running Tests

The commands below are examples; you should replace file paths and test names with your actual targets.

- **Run all tests in a file:**
    `pytest test/structure/test_structure.py`
- **Run a specific test function:**
    `pytest test/structure/test_structure.py::test_deepcopy`
- **Run tests by name pattern:**
    `pytest -k "deepcopy"`
- **Run tests by marker:**
    `pytest -m "interface"`

IMPORTANT! Always run the unit test you wrote and make sure they passes so that the function is complete and working!

#### Test Markers

The following markers are available to run specific subsets of tests:

- `accre`: Tests that should only run on the ACCRE cluster.
- `long`: Time-consuming tests (>10 minutes) that may need to be submitted as jobs.
- `interface`: Tests for the main software interfaces.
- `temp`: Temporary tests created during development.
- `TODO`: Tests that are known to be failing.

*Note: Not all eligible tests that fit a marker type are marked, but every marked test is eligible.*

## Architecture Overview

EnzyHTP is a modular Python library for high-throughput enzyme modeling that automates molecular simulation workflows. The architecture follows a layered design with clear separation of concerns.

### Core Architecture Principles
- **Modular Design**: Each workflow component is encapsulated in its own module
- **Interface Pattern**: Science APIs are decoupled from external tools via interface layer
- **Structure-Centric**: Most operations work on `Structure` or `StructureEnsemble` objects
- **Extensibility**: Easy to add new modules or replace existing implementations

### Key Modules

#### Core Data Structures
The library is built around several key data types that represent different aspects of molecular systems:

**Primary Structure Types (`enzy_htp/structure/`):**
- `Structure` - Central data type representing complete enzyme structures
- `StructureEnsemble` - Collections of structures for ensemble-based operations
- `StructureRegion` - Specific regions within structures for targeted operations

**Molecular Component Types:**
- `Atom`, `Residue`, `Chain` - Hierarchical molecular building blocks
- `Ligand`, `MetalAtom`, `Solvent` - Specialized molecular components
- `ModifiedResidue`, `NoncanonicalBase` - Non-standard structural elements

**Workflow Data Types:**
- `Mutation` (`enzy_htp/mutation_class/`) - Represents amino acid substitutions
- `ElectronicStructure` (`enzy_htp/electronic_structure/`) - Represents electronic structure and wavefunctions from quantum calculations
- Interface handle types in `enzy_htp/_interface/handle_types/`

**Structure Support Classes:**
- Structure I/O, operations, selections, constraints, and translation utilities
- Located in respective subdirectories under `enzy_htp/structure/`

Most functions in this library operate on these data types, particularly `Structure` objects.

#### Science APIs (High-level Workflow Modules)
These modules provide the main scientific functionality and are designed to be extensible:

**Current Science API Modules:**
- `enzy_htp/preparation/` - Structure preparation, cleaning, protonation, docking
- `enzy_htp/mutation/` - Amino acid substitutions and mutation pattern matching
- `enzy_htp/geometry/` - Reaction geometry and conformational sampling
- `enzy_htp/quantum/` - Quantum mechanics calculations and electronic structure
- `enzy_htp/analysis/` - Post-simulation analysis (RMSD, binding, stability, clustering, etc.)

**Workflow and Application Modules:**
- `enzy_htp/workflow_app/` - Complete workflow applications for specific use cases
- `enzy_htp/workflow_generation/` - Workflow configuration and execution management. (deporecated)

**Extensibility Note:**
New science API modules may be added in the future following the same patterns. Each typically contains:
- Submodules containing Science APIs (functions typically have "Science API" in their docstring)
- Clear separation from interface implementations

#### Interface Layer (`enzy_htp/_interface/`)
- Interfaces to external tools: Amber, Gaussian, Rosetta, PyMOL, RDKit, etc.
- Configurations stored in `enzy_htp/_config/`
- Science APIs use these interfaces but are not tightly coupled to them

#### Support Modules
- `enzy_htp/core/` - General utilities, file system, job management, clusters
- `enzy_htp/chemical/` - Knowledge base (chemistry, atoms, residues, physics constants)

## Code Conventions

### Code Quality

- MAKE SURE you write state-of-the-art quality code. Make it modular, robust, extensible, and maintainable.
- `yapf --in-place <file_path>`: Formats a Python file using the project's `yapf` style.
- **IMPORTANT:** Only use yapf for new files. For existing files, follow the established formatting to avoid unnecessary changes.

### Import Organization
Order imports by length within two sections:
```python
import os
import shutil
from typing import List
from subprocess import run

from .logger import _LOGGER
from .exception import MissingEnvironmentElement
```

### Key Patterns
- Functions primarily operate on `Structure` or `StructureEnsemble` objects
- Use interface pattern to avoid tight coupling with external tools
- Test-driven development with tests mirroring source code structure
- Type hints for function parameters and return values
- Google Python Style Guide compliance with custom pylint configuration
