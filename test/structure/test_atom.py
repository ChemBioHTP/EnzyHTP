import numpy as np
from enzy_htp.structure import Atom







def test_atom_distance():

	origin = Atom(x_coord=0,y_coord=0,z_coord=0)
	for (x_coord,y_coord,z_coord) in [
					[1,1,1],
					[-1,1,1],
					[1,-1,1],
					[1,1,-1],
					[-1,-1,1],
					[1,-1,-1],
					[-1,-1,-1],
					]:
		other = Atom(x_coord=x_coord,y_coord=y_coord,z_coord=z_coord)

		assert np.isclose(origin.distance( other ), np.sqrt(3) )


	origin = Atom(x_coord=0,y_coord=0,z_coord=0)
	other = Atom(x_coord=0,y_coord=0,z_coord=0)

	assert np.isclose(origin.distance( other ), 0 )


def test_atom_substract():
	pass
