import pytest


class Helpers:
    """helper functions for unit tests in EnzyHTP"""

    @staticmethod
    def lines_from_file(fname):
        fh = open(fname, "r")
        result = fh.read().splitlines()
        fh.close()
        return result

    @classmethod
    def equiv_files(cls, fname1: str, fname2: str, width: int = None,
                    skip_frist: bool = False, consider_order: bool= True) -> bool:
        """Helper method to check if two files are exactly equivalent."""
        if consider_order:
            first = 1
            for l1, l2 in zip(cls.lines_from_file(fname1), cls.lines_from_file(fname2)):
                if skip_frist and first:
                    first = 0
                    continue
                if width:
                    l1 = l1[:width]
                    l2 = l2[:width]

                if l1 != l2:
                    print(f"'{l1}'")
                    print(f"'{l2}'")
                    return False
            return True
        else:
            l1 = set(cls.lines_from_file(fname1))
            l2 = set(cls.lines_from_file(fname2))
            return l1 == l2


@pytest.fixture
def helpers():
    return Helpers
