from __future__ import division
import numpy as np
import time

import sidm_orbit_calculation.src.potentials.density as density
from sidm_orbit_calculation.src.utils.constants import *

class GetDragForce :
    def __init__(self, host) :
        self.host = host

        self.density_function = density.density_dict[host.potential]

        self.force = [0., 0., 0.]
        self.force_array = []

    def calculate_drag_force(self, position, momentum) :
        sigma_mDM = 3 * (1e-2*m_to_Mpc)**2 * M_sol # cm**2/g --> Mpc**2/M_sol # number from kalhoefer paper
        # sigma_mDM = 10 # testing
        # self.host.update_density(position)
        self.host.rho = self.calculate_density(position)
        vec  = 0.25 * sigma_mDM * momentum**2 * self.host.rho
        vec *= -np.sign(momentum)
        
        '''
        print 'drag force calculation'
        print 'pos = ', position
        print 'mom = ', momentum
        print 'force = ', vec
        print 
        time.sleep(2)
        '''

        return vec

    def calculate_density(self, position):
        return self.density_function(position[0],position[1],position[2],self.host)