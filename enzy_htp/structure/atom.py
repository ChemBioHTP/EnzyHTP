import numpy as np


class Atom:
	
	def __init__(self, **kwargs ):
		self.name = kwargs.get('name',str())
		self.x = kwargs.get('x',float())
		self.y = kwargs.get('y',float())
		self.z = kwargs.get('z',float())

	
	def __add__(self, offset):
		return Atom(
			name=self.name,
			x=self.x+other,
			y=self.y+other,
			z=self.z+other
			)

	def __sub__(self, offset):
		return Atom(
			name=self.name,
			x=self.x-other,
			y=self.y-other,
			z=self.z-other
			)


	def distance(self, other):
		return np.sqrt(
			(self.x-other.x)**2 +
			(self.y-other.y)**2 +
			(self.z-other.z)**2 
		)
