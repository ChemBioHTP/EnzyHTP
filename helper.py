'''
Misc helper func and class
'''
import os
import numpy as np
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
        os.makedirs(dir)


'''
geom
'''
def set_distance(p1,p2,d):
    '''
    return a coord of p3
    -- p3 --
    origin:     p1
    direction: (p1,p2)
    distance:   d
    '''
    p1 = np.array(p1)
    p2 = np.array(p2)
    #direction vector
    v1 = (p2 - p1)/np.linalg.norm(p1 - p2)
    p3 = p1 + v1 * d

    return tuple(p3)