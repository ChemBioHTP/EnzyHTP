# Notes about development


## Helper Scripts

We have developed a few helper scripts to aid developers when contributing to EnzyHTP. 
All scripts should be run from the root directory that `enzy_htp` is installed in.
Their descriptions
are seen below:

+ `install`: installs `enzy_htp` from source using `pip`. Detects `pip` version and uses
appropriate commandline arguments.

+ `lint`: Runs the `yapf` formatter with the custom EnzyHTP settings on both the source code
found in `enzy_htp/` as well as the tests found in `test/`

+ `PR-prep`: Should be run before all pull request (PR's) are made. This script installs `enzy_htp` from source, 
runs unit tests and checks that the package can be installed from source correctly.




## python imports

When importing at the top of python files, put the non-`enzy_htp` imports at the top of and then the
`enzy_htp` imports below, with one line separating them. Make sure that the lines are ordered by increasing
length for **both** subsections. See example below:

```
import os
import shutil
from typing import List
from subprocess import run

from .logger import _LOGGER
from .exception import MissingEnvironmentElement

```
