'''
Misc helper func and class
'''
import os
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

'''
func
'''
def mkdir(dir):
    if os.path.exists(dir):
        pass
    else:
        os.mkdir(dir)
