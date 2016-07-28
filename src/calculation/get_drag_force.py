from __future__ import division
import numpy as np

import sidm_orbit_calculation.src.potentials.density as density

class GetDragForce :
    def __init__(self, host) :
        self.host = host

        self.density_function = density.density_dict[host.density]

        self.force = [0., 0., 0.]
        self.force_array = []

    def calculate_drag_force(self, position, momentum) :
        sigma_mDM = 1.e-4 # cm2/g 
        # sigma_mDM = 1e4 # testing
        # self.host.update_density(position)
        self.calculate_density(position)
        return 0.25 * sigma_mDM * momentum**2 * self.host.rho

	def calculate_density(self, position):
		self.host.rho = self.density_function(position[0],position[1],position[2],self.host)