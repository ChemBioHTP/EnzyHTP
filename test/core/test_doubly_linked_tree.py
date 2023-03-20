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
