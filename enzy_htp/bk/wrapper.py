"""
wappers and context manager defined for the program
"""
import os, sys

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
