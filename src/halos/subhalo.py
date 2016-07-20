import numpy as np

from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.calculation.get_drag_force import *
from sidm_orbit_calculation.src.orbit_parameters.orbit_distributions import *
from sidm_orbit_calculation.src.utils.constants import *

class Subhalo:

	def __init__(self, host, M, initial_position, initial_momentum):
		self.host = host
		self.M = M*M_sol
		self.initial_parameters(initial_position,initial_momentum)
		
		# forces should not necessarily be classes
		self.gravity = GetGravitationalForce(self.host)
		self.drag = GetDragForce(self.host)

		self.count = 0

	def initial_parameters(self,initial_position,initial_momentum):
		if not initial_position.any():
			self.position, self.momentum = initial_conditions(self)
		else:
			self.position = initial_position
			self.momentum = -initial_momentum
		print 'initial position = %.2e, %.2e, %.2e' % (self.position[0], self.position[1], self.position[2])
		print 'initial momentum = %.2e, %.2e, %.2e' % (-self.momentum[0], -self.momentum[1], -self.momentum[2])

		self.position_array = [(self.position[0],self.position[1],self.position[2])]
		self.momentum_array = [(-self.momentum[0],-self.momentum[1],-self.momentum[2])]


       
