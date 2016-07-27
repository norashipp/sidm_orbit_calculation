import numpy as np
from scipy.misc import derivative
import time

import sidm_orbit_calculation.src.potentials.test_spherical_potentials as potentials

'''This will need the gravitational potential to calculate the force
on a test particle at any point'''

class GetGravitationalForce:
    def __init__(self, host):
        self.potential = potentials.potential_dict[host.potential]
        # self.potential_function = lambda x,y,z: self.potential(x=x,y=y,z=z,host=host)
        self.potential_function = lambda r: self.potential(r=r,host=host)
        self.force_array = []
        
    def calculate_gravitational_force(self, position):
        partial = 0
        if partial:
            mag, vec = self.calculate_partial_force(position)
        else:
            mag, vec = self.calculate_force(position)
        return mag, vec

    def calculate_force(self,position):
        dx = 1e20 # is this causing problems?
        r, theta, phi = spherical_coordinates(position)
        mag = derivative(self.potential_function,r,dx)
        vec = self.get_vector(mag,phi,theta)
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
        return mag, vec

    def get_vector(self, mag, phi, theta):
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

def spherical_coordinates(position):
    r = np.sqrt(position[0]**2+position[1]**2+position[2]**2)
    theta = np.arccos(position[2]/r)
    phi = np.arctan2(position[1],position[0])
    return r, theta, phi