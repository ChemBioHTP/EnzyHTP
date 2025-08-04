# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

### Installation and Setup
- Avoid installing the package yourself. If you found a problem involving environment or need to test install, ask for help from the user.
- `./dev-tools/pip-install` - Install enzy_htp from source using pip (do not install dependencies)
- `./dev-tools/conda-install` - Install enzy_htp and dependencies using conda

### Code Quality and Testing
- `./dev-tools/lint` - Run yapf formatter on source code and tests (uses custom EnzyHTP settings) NOTE! only use this in new files.
- `pytest` - Run unit tests (supports markers: `accre`, `long`, `clean`, `interface`, `temp`, `TODO`)
- `pytest -m "not long"` - Run tests excluding time-consuming ones
- `pytest -m accre` - Run tests that require ACCRE cluster environment

### Test Markers
- `accre`: Tests that only run on ACCRE cluster
- `long`: Time-consuming tests (>10min) that need job submission
- `interface`: Tests for main interfaces
- `temp`: Temporary development tests
- `TODO`: Known failing tests
NOTE Not all eligible tests that fit a maker type are marked, but every marked test is eligible.

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

### Code Conventions

#### Import Organization
Order imports by length within two sections:
```python
import os
import shutil
from typing import List
from subprocess import run

from .logger import _LOGGER
from .exception import MissingEnvironmentElement
```

#### Key Patterns
- Functions primarily operate on `Structure` or `StructureEnsemble` objects
- Use interface pattern to avoid tight coupling with external tools
- Test-driven development with tests mirroring source code structure
- Type hints for function parameters and return values
- Google Python Style Guide compliance with custom pylint configuration