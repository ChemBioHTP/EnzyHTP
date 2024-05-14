from __future__ import annotations
from typing import Tuple
from dataclasses import dataclass

from .residue import convert_to_one_letter, THREE_LETTER_AA_MAPPER

@dataclass
class SeqRes:
    
    model:int
    chain:str
    idx:int
    name:str
    missing:bool
    seq_idx:int

    def one_letter(self) -> str:

        return convert_to_one_letter(self.name)

    @property
    def key(self) -> Tuple[str, int]:
        return (self.chain, self.idx)

    def is_canonical(self) -> bool:
        return self.name in THREE_LETTER_AA_MAPPER
        

    def __gt__(self, other) -> bool:
        return self.key > other.key


    def __lt__(self, other) -> bool:
        return self.key < other.key


    
