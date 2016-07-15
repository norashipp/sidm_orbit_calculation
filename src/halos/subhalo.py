import numpy as np

from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.calculation.get_drag_force import *
from sidm_orbit_calculation.src.orbit_parameters.orbit_distributions import *

class Subhalo:

	def __init__(self, host, M, initial_position, initial_momentum):
		self.host = host
		self.M = M
		self.initial_parameters(initial_position,initial_momentum)
		
		# forces should not necessarily be classes
		self.gravity = GetGravitationalForce(self.host)
		self.drag = GetDragForce(self.host)

	def initial_parameters(self,position,momentum):
		if position and momentum:
			self.position = initial_position
			self.momentum = initial_momentum
		else:
			self.position, self.momentum = initial_conditions(self)
			print 'initial position = ', self.position
			print 'initial momentum = ', self.momentum

		self.position_array = [(self.position[0],self.position[1],self.position[2])]
		self.momentum_array = [(self.momentum[0],self.momentum[1],self.momentum[2])]


       
