'''
Misc helper func and class
'''

'''
====
Tree
====
'''
class Child():
    def __init__(self):
        self.parent = None
    def set_parent(self, parent_obj):
        self.parent = parent_obj

        return self

'''
Text
'''
line_feed = '\n'