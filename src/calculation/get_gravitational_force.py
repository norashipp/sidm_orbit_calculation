import numpy as np
from scipy.misc import derivative
import time

import sidm_orbit_calculation.src.potentials.test_spherical_potentials as potentials

'''This will need the gravitational potential to calculate the force
on a test particle at any point'''

class GetGravitationalForce :
    def __init__(self, host) :
        self.potential = potentials.potential_dict[host.potential]
        self.potential_function = lambda r: self.potential(r=r,host=host)
        self.force_array = []

    def calculate_gravitational_force(self, position) :
        r = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
        dx = 1.e20
        mag = derivative(self.potential_function,r,dx)
        return mag

    def get_vector(self, mag, phi,theta):
        fx = -1.*mag*np.cos(phi)*np.sin(theta)
        fy = -1.*mag*np.sin(phi)*np.sin(theta)
        fz = -1.*mag*np.cos(theta)
        vec = np.array([fx, fy, fz])
        print 'gravity'
        print vec
        print
        time.sleep(2)
        return vec