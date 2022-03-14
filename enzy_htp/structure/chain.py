# TODO what the heck else does this class do???
import logging
from typing import List

from .residue import Residue

class Chain:
    def __init__(self, name, residues):
        self.name_ = name
        self.residues_ : List[Residue] = residues
        # CJ: not sure if I should do this but I am overwriting every time
        # right now
        _ = list(map(lambda r: r.set_chain(self.name_), self.residues_))

    def is_metal(self) -> bool:
        for rr in self.residues_:
            if rr.is_metal():
                logging.warn(f"\033[1;34;40mStructure: found metal in raw: {self.name_} {rr.name} {residue.id} \033[0m")
                return True
        return False
