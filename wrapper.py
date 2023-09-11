'''
wappers and context manager defined for the program
'''
import os,sys
import logging


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
        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr
        self.redirect_stdout = open(self.redirect, 'w')
        self.redirect_stderr = open(self.redirect, 'w')
        sys.stdout = self.redirect_stdout
        sys.stderr = self.redirect_stderr

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.redirect_stdout.close()
        self.redirect_stderr.close()
        sys.stdout = self.original_stdout
        sys.stderr = self.original_stderr
        # restore any logging handler that lead to redirect_io_obj
        all_loggers = [logging.getLogger("")] + [j for i, j in logging.Logger.manager.loggerDict.items()]
        for logger in all_loggers:
            if hasattr(logger, "handlers"):
                for handler in logger.handlers:
                    if hasattr(handler, "stream"):
                        if handler.stream is self.redirect_stdout:
                            handler.stream = self.original_stdout
                        if handler.stream is self.redirect_stderr:
                            handler.stream = self.original_stderr
