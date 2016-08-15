import numpy as np

from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.calculation.get_drag_force import *
from sidm_orbit_calculation.src.orbit_parameters.orbit_distributions import *
from sidm_orbit_calculation.src.utils.constants import *

class Subhalo:

	def __init__(self, host, M, initial_position, initial_momentum, t0):
		self.host = host
		self.t0 = t0
		self.M = M # *M_sol # work in solar masses
		self.initial_parameters(initial_position,initial_momentum)
		
		# forces should not necessarily be classes
		self.gravity = GetGravitationalForce(self.host)
		self.drag = GetDragForce(self.host)

		self.count = 0

		self.energy_array = []
		
	def initial_parameters(self,initial_position,initial_momentum):
		if not initial_position.any():
			self.position, self.momentum = initial_conditions(self)
			self.position *= self.host.R # initial position in units of host halo virial radius
		else:
			self.position = initial_position # *self.host.R
			self.momentum = initial_momentum # *self.host.v # *m_to_Mpc/s_to_Gyr

		# print 'initial position = %.2e, %.2e, %.2e' % (self.position[0], self.position[1], self.position[2])
		# print 'initial momentum = %.2e, %.2e, %.2e' % (self.momentum[0], self.momentum[1], self.momentum[2])

		self.position_array = [(self.position[0],self.position[1],self.position[2])]
		self.momentum_array = [(self.momentum[0],self.momentum[1],self.momentum[2])]

	def calculate_energy(self):
		KE = 0.5*self.M*(self.momentum[0]**2+self.momentum[1]**2+self.momentum[2]**2)
		PE = self.gravity.potential_function(self.position[0], self.position[1], self.position[2])
		E = KE + PE
		self.energy_array.append(E)