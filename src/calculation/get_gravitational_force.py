import numpy as np
from scipy.misc import derivative
import sidm_orbit_calculation.src.potentials.test_spherical_potentials as potentials

'''This will need the gravitational potential to calculate the force
on a test particle at any point'''

class GetGravitationalForce :
    def __init__(self, host=None) :
        self.potential = potentials.potential_dict[host.potential]
        self.potential_function = lambda r: self.potential(r=r,host=host)

        self.force_array = []

    def calculate_gravitational_force(self, position=None) :
        dx = position[0]/42 # gives smallest final r
        # force_r = derivative(self.potential_function,position[0],dx)
        force_r = partial_derivative(self.potential_function,0,[position[0]], dx=dx) # generalize (!!)
        force_phi = 0. # for now just spherical potentials (!!)
        self.force_array.append(force_r)
        return force_r,force_phi

    def partial_force(self,position=None):
        dx = position[0]/42
        partial_force_r = derivative(self.potential_function,position[0],dx,n=2) # works as long as potential is spherically symmetric
        # force = lambda r: partial_derivative(self.potential_function,0,r, dx=dx)
        # partial_force_r = partial_derivative(force,0,[position[0]], dx=dx)
        return partial_force_r
#################################################################    

def partial_derivative(func, var=0, point=[], dx=1.e10):
        args =  point[:]
        def wraps(x):
            args[var] = x
            return func(*args)
        return derivative(wraps, point[var], dx = dx)