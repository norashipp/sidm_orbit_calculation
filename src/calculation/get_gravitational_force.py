import numpy as np
from scipy.misc import derivative
import time

import sidm_orbit_calculation.src.potentials.test_spherical_potentials as potentials

'''This will need the gravitational potential to calculate the force
on a test particle at any point'''

class GetGravitationalForce:
    def __init__(self, host):
        self.potential = potentials.potential_dict[host.potential]
        self.potential_function = lambda x,y,z: self.potential(x=x,y=y,z=z,host=host)
        self.force_array = []
        
        # probably a better way to do this... (!!)
        # spherical = 1
        # if host.potential == 'triaxial_NFW':
        #     spherical = 0

    def calculate_gravitational_force(self, position):
        # if not spherical:
        mag = self.calculate_partial_force(position)
        # else:
        #     mag = derivative(self.potential_function,r,dx)
        return mag

    def calculate_partial_force(self,position):
        dx = 1e16
        dpx = partial_derivative(self.potential_function,0,position,dx)
        dpy = partial_derivative(self.potential_function,1,position,dx)
        dpz = partial_derivative(self.potential_function,2,position,dx)
        return dpx + dpy + dpz

    def get_vector(self, mag, phi,theta):
        fx = -1.*mag*np.cos(phi)*np.sin(theta)
        fy = -1.*mag*np.sin(phi)*np.sin(theta)
        fz = -1.*mag*np.cos(theta)
        vec = np.array([fx, fy, fz])
        ep = 1e-16
        vec[np.where(np.abs(vec) <= 1000*ep)] = 0

        # print vec
        # print
        # time.sleep(2)
        
        return vec

def partial_derivative(func, var=0, point=[], dx=1.e10):
    args =  point[:]
    def wraps(x):
        args[var] = x
        return func(*args)
    return derivative(wraps, point[var], dx = dx)