import numpy as np
from scipy.misc import derivative
from scipy.integrate import quad
import time

import sidm_orbit_calculation.src.potentials.test_spherical_potentials as potentials
from sidm_orbit_calculation.src.utils.geometry import *

'''This will need the gravitational potential to calculate the force
on a test particle at any point'''

class GetGravitationalForce:
    def __init__(self, host):
        self.potential = potentials.potential_dict[host.potential]
        self.potential_function = lambda x,y,z: self.potential(x=x,y=y,z=z,host=host)
        
        self.force_array = []
        # self.count = 0
        # self.t0 = time.time()

    def calculate_density(self,position):
        return self.density_function(position[0], position[1], position[2], host)

    def calculate_gravitational_force(self, position):
        mag, vec = self.calculate_partial_force(position)
        return mag, vec
           
    def calculate_partial_force(self, position):
        dx = 1e18 # is this causing problems?
        fx = partial_derivative(self.potential_function,0,position,dx)
        fy = partial_derivative(self.potential_function,1,position,dx)
        fz = partial_derivative(self.potential_function,2,position,dx)
        
        vec = np.array([-fx,-fy,-fz])
        ep = 1e-16
        vec[np.where(np.abs(vec) <= 1000*ep)] = 0
        
        mag = np.sqrt(fx**2 + fy**2 + fz**2)

        # self.count+=1
        # if self.count % 100 == 0:
        #     t1 = time.time()
        #     print self.count, t1 - self.t0
        #     self.t0 = t1
        
        return mag, vec

def partial_derivative(func, var=0, point=[], dx=1.e10):
    args =  point[:]
    def wraps(x):
        args[var] = x
        return func(*args)
    return derivative(wraps, point[var], dx = dx)
