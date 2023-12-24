"""This class summarized data structure level operation from Structure/Chain/Residue/Atom, which
are designed to contain data in the doublelinked manner, that is, a parent object holds a list
of children objects while each children object also holds a reference of the parent object. (e.g:
Chain().residues[0].chain is the first Chain object)

Common operations from such summarization are:
set/get_parent (set_ghost_children for leafs)
set/get_children (set_ghost_parent for the root, these are designed so that the summarized class
                  know what to expect.)
__deepcopy__
__getitem__
__delitem__
__len__

Author: QZ Shao, <shaoqz@icloud.com>
Date: 2022-09-14
"""

import copy
from typing import Any, Dict, List, Union
from .general import delete_base_on_id


class DoubleLinkedNode():
    """
    class for parent objects of the doubly linked tree
    """

    def __init__(self, children: List = None, parent=None):
        """place holder init function for deepcopy testing"""
        if children is None:
            self.set_ghost_children()
        else:
            self.set_children(children)
        if parent is None:
            self.set_ghost_parent()
        else:
            self.set_parent(parent)

    #region === Attr ===
    # parent use
    def set_children(self, children: List):
        """
        set children and add self as parent of children
        """
        self._children = children
        for child in self._children:
            child.parent = self

    def set_ghost_children(self):
        """
        method for the node with no children to set
        The idea is DoubleLinkNode defines what to set when no children
        """
        self._children = None

    def get_children(self) -> List:
        return self._children

    # api out of class use
    @property
    def children(self) -> List:
        return self.get_children()

    @children.setter
    def children(self, val):
        self.set_children(val)

    # child use
    def set_parent(self, parent):
        """
        only set parent to self
        * do not add self to parent"s children
        """
        self._parent = parent

    def set_ghost_parent(self):
        """
        method for the node with no parent to set
        The idea is DoubleLinkNode defines what to set when no parent
        """
        self._parent = None

    def get_parent(self):
        return self._parent

    # api out of class use
    @property
    def parent(self):
        return self.get_parent()

    @parent.setter
    def parent(self, val):
        self.set_parent(val)

    @property
    def root(self):
        """get the very parent object that has no parent"""
        if self.parent is None:
            return self
        return self.parent.root

    #endregion

    #region === edit ===
    def delete_from_parent(self) -> None:
        """
        delete self from parent's children list
        delete base on object's id()
        Returns:
            None. changes are made to the .children list in self.parent
        """
        delete_base_on_id(self.parent.children, id(self))

    def __deepcopy__(self, memo: Union[Dict[int, Any], None] = None, _nil=[]):
        """
        Support deepcopy of DoublyLinkedNode that donot copy any parent and siblings.
        (Implemtation inspired by https://stackoverflow.com/a/40484215)
        The default deepcopy method copies every object in the tree. (e.g.: copying 3
        atoms one by one will generate 3 new complete structure, which will be much slower
        and redundant.) Following the suggestion of the doc string of the copy module,
        here re-define the deepcopy behavior by not copying the parent object and set it
        to None. (copying of the entrie tree should be done by copying the root.)

        Returns:
            a deepcopy of a DoubleLinkedNode with parent = None.

        NOTE: if there are not constructor in the class. All children with have no parent.
        """
        # in case this is the first copied item
        if memo is None:
            memo = {}
        # treat copying action on parent
        if self.parent is not None:
            # not the root.
            parent_id = id(self.parent)
            #print(f"parent_id: {parent_id}")
            y = memo.get(parent_id, _nil)
            if y is _nil:
                # parent not in memo -> this is the first copied DoublyLinkNode: set parent to None in memo
                memo[id(self.parent)] = None
            # parent in memo -> this is part of the recursive copying initiated from the parent: default

        # mask current method to use original deepcopy
        self.__deepcopy__ = None

        new_self = copy.deepcopy(self, memo)

        # recover current method for future use
        # test shows this will rebind the method to correct instance when its called again
        delattr(self, "__deepcopy__")
        delattr(new_self, "__deepcopy__")

        return new_self

    # def deepcopy_complete_tree(self, memo=None, _nil=[]):
    #     """this reproduct the default deepcopy behavior. Not working though. Since it will not be called recursively"""
    #     # Start memo from the first item
    #     if memo is None:
    #         memo = {}
    #     # give reference of the correponding new self when a self is already copied once.
    #     d = id(self)
    #     y = memo.get(d, _nil)
    #     if y is not _nil:
    #         return y
    #     # actually copy self if never copied before
    #     reductor = getattr(self, "__reduce_ex__", None)
    #     func, args, listiter, dictiter = reductor(4)
    #     # _reconstruct
    #     if args:
    #         args = (copy.deepcopy(arg, memo) for arg in args)
    #     y = func(*args)
    #     memo[id(self)] = y

    #     if state is not None:
    #         state = copy.deepcopy(state, memo)
    #         if isinstance(state, tuple) and len(state) == 2:
    #             state, slotstate = state
    #         else:
    #             slotstate = None
    #         if state is not None:
    #             y.__dict__.update(state)
    #         if slotstate is not None:
    #             for key, value in slotstate.items():
    #                 setattr(y, key, value)
    #     memo[d] = y
    #     # _keep_alive
    #     try:
    #         memo[id(memo)].append(self)
    #     except KeyError:
    #         # aha, this is the first one :-)
    #         memo[id(memo)]=[self]
    #     return y
    #endregion

    # === special ===
    def __getitem__(self, key: int):
        return self._children[key]

    def __delitem__(self, key: int):
        del self._children[key]

    def __len__(self) -> int:
        return len(self._children)
