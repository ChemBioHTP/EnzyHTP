'''
wappers and context manager defined for the program
'''
import os,sys


# def blockprint(redirect=os.devnull):
#     '''
#     blockprint(redirect=os.devnull)
#     -----
#     block or redirect stdout prints to "redirect"
#     '''
#     def blockprint_dec(func):
#         def wrapper(*args, **kwargs):
#             sys.stdout = open(redirect, 'w')
#             results = func(*args, **kwargs)
#             sys.stdout = sys.__stdout__
#             return results
#         return wrapper
#     return blockprint_dec

class HiddenPrints:
    '''
    HiddenPrints(redirect=os.devnull)
    -----
    block or redirect stdout prints to "redirect"
    '''
    def __init__(self, redirect=os.devnull):
        self.redirect=redirect

    def __enter__(self):
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = open(self.redirect, 'w')
        sys.stderr = open(self.redirect, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr
