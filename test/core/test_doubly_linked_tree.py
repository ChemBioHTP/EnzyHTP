"""Testing the enzy_htp.core.DoubleLinkedNode() object.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-22
"""

import copy
import pytest

from enzy_htp.core.doubly_linked_tree import DoubleLinkedNode
from enzy_htp.core.general import EnablePropagate
from enzy_htp.core.logger import _LOGGER


def test_delete_from_parent():
    # model objects
    child_1 = DoubleLinkedNode()
    child_2 = DoubleLinkedNode()
    parent = DoubleLinkedNode(children=[child_1, child_2])

    child_1.delete_from_parent()

    assert id(parent.children[0]) == id(child_2)


def test_deepcopy():
    """test if deepcopy assign parent correctly"""
    child_1 = DoubleLinkedNode()
    child_2 = DoubleLinkedNode()
    parent = DoubleLinkedNode(children=[child_1, child_2])

    new_parent = copy.deepcopy(parent)
    assert id(new_parent.children[0].parent) == id(new_parent)

def test_deepcopy_more_than_once():
    """test if deepcopy act as expected when it is copied more than once"""
    child_1 = DoubleLinkedNode()
    child_2 = DoubleLinkedNode()
    parent = DoubleLinkedNode(children=[child_1, child_2])
    parent[0].mark_temp_1 = 1
    assert hasattr(parent[0], "mark_temp_1")
    assert not hasattr(parent[0], "mark_temp_2")

    # 1st copy
    new_parent = copy.deepcopy(parent)
    new_parent[0].mark_temp_2 = 1
    assert id(new_parent) != id(parent)
    assert hasattr(new_parent[0], "mark_temp_1")
    assert hasattr(new_parent[0], "mark_temp_2")

    # 2nd copy (make sure it is not copying the original object)
    new_new_parent = copy.deepcopy(new_parent)
    assert id(new_new_parent) != id(parent)
    assert id(new_new_parent) != id(new_parent)
    assert hasattr(new_new_parent[0], "mark_temp_1")
    assert hasattr(new_new_parent[0], "mark_temp_2")

def test_root():
    """test using a made up tree."""
    child_1 = DoubleLinkedNode()
    child_2 = DoubleLinkedNode()
    parent_1 = DoubleLinkedNode(children=[child_1, child_2])
    child_3 = DoubleLinkedNode()
    child_4 = DoubleLinkedNode()
    parent_2 = DoubleLinkedNode(children=[child_3, child_4])
    root = DoubleLinkedNode(children=[parent_1, parent_2])

    assert child_1.root() is root

def test_root_inf_loop(caplog):
    """test using a made up tree containing a inf loop"""
    child_1 = DoubleLinkedNode()
    child_2 = DoubleLinkedNode()
    parent_1 = DoubleLinkedNode(children=[child_1, child_2])
    child_3 = DoubleLinkedNode()
    child_4 = DoubleLinkedNode()
    parent_2 = DoubleLinkedNode(children=[child_3, child_4])
    root = DoubleLinkedNode(children=[parent_1, parent_2])
    root.parent = child_1

    with EnablePropagate(_LOGGER):
        with pytest.raises(RuntimeError) as e:
            child_1.root()
            assert "inf" in caplog.text

