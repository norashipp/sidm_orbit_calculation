import numpy as np
from scipy.misc import derivative
import sidm_orbit_calculation.src.potentials.test_spherical_potentials as potentials

'''This will need the gravitational potential to calculate the force
on a test particle at any point'''

class GetGravitationalForce :
    def __init__(self, host=None) :
        self.potential = potentials.potential_dict[host.potential]
        self.potential_function = lambda r: self.potential(r=r,host=host)
        self.force = 0.
        self.force_array = []

    def calculate_gravitational_force(self, position=None) :
        r = np.sqrt(position[0]**2+position[1]**2)
        dx = 1.e20
        # dx = position[0]/42 # gives smallest final r
        self.force = derivative(self.potential_function,r,dx)
        # force = partial_derivative(self.potential_function,0,r, dx=dx) # generalize (!!)
        self.force_array.append(self.force)
 
    def calculate_partial_force(self,position=None):
        dx = 1.e20
        partial_force = derivative(self.potential_function,position[0],dx,n=2) # works as long as potential is spherically symmetric
        return np.array([partial_force,0]) # generalize (!!)

    def vector(self,phi=None):
        return np.array([-1.*self.force*np.cos(phi), -1.*self.force*np.sin(phi)])

def partial_derivative(func, var=0, point=[], dx=1.e10):
        args =  point[:]
        def wraps(x):
            args[var] = x
            return func(*args)
        return derivative(wraps, point[var], dx = dx)