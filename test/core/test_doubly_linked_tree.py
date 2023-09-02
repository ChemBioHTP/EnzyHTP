"""Testing the enzy_htp.core.DoubleLinkedNode() object.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-22
"""

import copy

from enzy_htp.core.doubly_linked_tree import DoubleLinkedNode


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
