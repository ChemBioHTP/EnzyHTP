'''
helper for check whats a pandas df looks like
'''

from biopandas.pdb import PandasPdb
import sys

test_mdl = sys.argv[1]
test_mdl_pdb = PandasPdb()
test_mdl_pdb.read_pdb(test_mdl)
print(test_mdl_pdb.df['ATOM'])
