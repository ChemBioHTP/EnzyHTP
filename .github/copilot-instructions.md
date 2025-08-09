# Project Overview

This project is a Python library that automates the building blocks of computational enzyme modeling workflows by converting each modular step into Python APIs. This is achieved by interfacing with popular command-line software as well as developing new algorithms.

## Folder Structure

- `/enzy_htp`: Contains the source code of this project.
- `/test`: Contains unit tests of the source code based on the pytest framework.
- `/template`: Contains template files of example workflows used in the tutorial.
- `/dev-tools`: Contains handy tool scripts used during install or development
- `/resources`: Contains references or figures as supporting information.
- `/test_integration`: Contains integration tests of example workflows. (only work in older version; deprecated)

## Architecture

The architecture of this project is designed to be modular and extensible. Each component of the enzyme modeling workflow is encapsulated within its own module, allowing for easy integration and customization. The core data structures are defined in the `enzy_htp/structure` module, while the high-level Science APIs are implemented in their respective modules under `enzy_htp/preparation`, `enzy_htp/mutation`, etc. These Science APIs are supported, but not coupled, with interfaces to popular command-line software or Python packages, which are defined in the `enzy_htp/_interface` module. The configurations for these interfaces are stored in the `enzy_htp/_config` module. The Science APIs are not coupled with the interfaces, allowing users to choose their preferred "engine" for each step of the workflow. This modular design enables users to easily extend the library by adding new modules or replacing existing ones with custom implementations.

Highlighted modules in the architecture:
- `enzy_htp/structure` contains core data structures. The most important one is `Structure`, which represents an enzyme structure. Most functions in this library operate on `Structure` objects. (other important ones are `StructureEnsemble`, `StructureRegion`, etc.)
- `enzy_htp/_interface` contains interfaces to popular command-line software or python packages.
- `enzy_htp/_config` contains configurations of each interfaced software or package.
- `enzy_htp/preparation`, `enzy_htp/mutation`, `enzy_htp/geometry`, `enzy_htp/quantum`, `enzy_htp/structure_prediction`, and `enzy_htp/analysis` contain high-level APIs of modular steps of the enzyme modeling workflow. We call them "Science APIs".
- `enzy_htp/core` contains general utilities and helper functions that are used across the library.
- `enzy_htp/chemical` contains knowledge about chemistry, physics, etc. (will rename to `enzy_htp/knowledge` in the future)

## Coding Standards

- Always think about future extensibility and modularity when writing code.
- Write unit tests for complicated functions.
- Always use existing APIs or libraries to implement a function, rather than writing it from scratch.
- Follow Python style conventions (PEP 8) and Google Python Style Guide.
- Use type hints for function parameters and return values.
- Handle exceptions appropriately and provide meaningful error messages.

## Key Patterns

- Most functions operate on `Structure` object or `StructureEnsemble` object as the primary data type.
- Use the interface pattern: Science APIs should not be tightly coupled with specific external tools.
- Test-driven development: Write tests in the `/test` directory following the same structure as source code.
