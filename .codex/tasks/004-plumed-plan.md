# PLUMED-AMBER Integration Implementation Plan

## Overview
This plan outlines the steps required to support running PLUMED with AMBER in EnzyHTP. PLUMED is a plugin for enhanced sampling in molecular dynamics, and AMBER is a widely used MD engine. The goal is to enable users to set up, run, and analyze AMBER simulations with PLUMED biasing via EnzyHTP's modular workflow.

## 1. Requirements
- **PLUMED**: Installed and compiled with AMBER support (patched AMBER).
- **AMBER**: AMBER20+ (pmemd, sander) with PLUMED patch applied.
- **EnzyHTP**: Modular interface to configure, launch, and analyze PLUMED-enabled AMBER jobs.

## 2. High-Level Steps
1. **PLUMED Input Preparation**
    - Generate PLUMED input files (e.g., `plumed.dat`) from user specifications or templates.
    - Allow user to specify collective variables (CVs), biasing methods, and output options.
2. **AMBER Input Preparation**
    - Ensure AMBER input files are compatible with PLUMED (e.g., correct &cntrl flags).
    - Add PLUMED-specific flags to AMBER input (e.g., `&cntrl ifplumed=1` for sander).
3. **Job Submission**
    - Modify job submission scripts to include PLUMED environment variables if needed.
    - Ensure correct AMBER executable is used (PLUMED-patched).
4. **Execution**
    - Run AMBER with PLUMED, ensuring `plumed.dat` is in the working directory.
    - Capture and log PLUMED output files.
5. **Analysis**
    - Parse PLUMED output (e.g., COLVAR, HILLS files) for downstream analysis.
    - Integrate with EnzyHTP analysis modules.

## 3. API/Interface Design
- **New/Updated Modules**:
    - `enzy_htp/_interface/plumed_interface.py`: Interface for PLUMED input generation and output parsing.
    - Update `enzy_htp/_interface/amber_interface.py` to support PLUMED flags and execution.
    - Update `enzy_htp/geometry` or `enzy_htp/workflow_app` for workflow integration.
- **Key Functions**:
    - `generate_plumed_input(structure, cv_spec, bias_spec, output_spec) -> str`
    - `run_amber_with_plumed(structure, amber_input, plumed_input, ...)`
    - `parse_plumed_output(output_files) -> dict`

## 4. User Workflow
1. User specifies CVs and biasing in EnzyHTP workflow (Python API or YAML/JSON config).
2. EnzyHTP generates `plumed.dat` and AMBER input files.
3. EnzyHTP launches AMBER with PLUMED.
4. EnzyHTP collects and parses results for analysis.

## 5. Testing
- Unit tests for PLUMED input generation and output parsing.
- Integration tests for running a minimal AMBER+PLUMED simulation (e.g., metadynamics on alanine dipeptide).

## 6. Documentation
- Usage examples in the EnzyHTP tutorial.
- Troubleshooting for common PLUMED-AMBER issues (e.g., patching, environment setup).

## 7. Future Extensions
- Support for other MD engines (GROMACS, NAMD) via the same PLUMED interface.
- Advanced CVs and biasing methods.

---
**References:**
- [PLUMED Documentation](https://www.plumed.org/doc-v2.8/user-doc/html/index.html)
- [AMBER Manual: PLUMED Support](https://ambermd.org/)
