# Notes about development


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
