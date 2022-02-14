import numpy as np
from enzy_htp.structure import Atom







def test_atom_distance():

	origin = Atom(x=0,y=0,z=0)
	for (x,y,z) in [
					[1,1,1],
					[-1,1,1],
					[1,-1,1],
					[1,1,-1],
					[-1,-1,1],
					[1,-1,-1],
					[-1,-1,-1],
					]:
		other = Atom(x=x,y=y,z=z)

		assert np.isclose(origin.distance( other ), np.sqrt(3) )


	origin = Atom(x=0,y=0,z=0)
	other = Atom(x=0,y=0,z=0)

	assert np.isclose(origin.distance( other ), 0 )


def test_atom_substract():
	pass
