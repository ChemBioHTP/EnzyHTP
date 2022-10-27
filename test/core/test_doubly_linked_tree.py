"""Testing the enzy_htp.core.DoubleLinkedNode() object.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2022-09-22
"""

from enzy_htp.core.doubly_linked_tree import DoubleLinkedNode


def test_delete_from_parent():
    # model objects
    child_1 = DoubleLinkedNode()
    child_2 = DoubleLinkedNode()
    parent = DoubleLinkedNode()

    parent.set_children([child_1, child_2])
    child_1.delete_from_parent()

    assert id(parent.children[0]) == id(child_2)
