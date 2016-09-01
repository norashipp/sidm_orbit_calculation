from __future__ import division

import numpy as np
from scipy.misc import derivative
from scipy.integrate import quad
import time
from numpy import linalg

import sidm_orbit_calculation.src.potentials.test_spherical_potentials as potentials
from sidm_orbit_calculation.src.utils.geometry import *

'''This will need the gravitational potential to calculate the force
on a test particle at any point'''

class GetGravitationalForce:
    def __init__(self, host):
        self.host = host
        self.potential = potentials.potential_dict[host.potential]
        self.potential_function = lambda x,y,z: self.potential(x=x,y=y,z=z,host=host)
        
        self.force_array = []
        # self.count = 0
        # self.t0 = time.time()

    def calculate_density(self,position):
        return self.density_function(position[0], position[1], position[2], host)

    def calculate_gravitational_force(self, position):
        position_rotated = self.rotate(position)
        mag, vec = self.calculate_partial_force(position_rotated)
        # mag, vec = self.calculate_partial_force(position)
        return mag, vec

    def rotate(self,position):
        # THIS ONLY WORKS WITH THE SINGLE AXIS DIRECTION
        axis = np.array([self.host.ax,self.host.ay,self.host.az])

        _, theta1, phi1 = spherical_coordinates(np.array([1,0,0]))
        _, theta2, phi2 = spherical_coordinates(axis)

        phi = np.pi/2 - phi1
        theta = theta1 - theta2
        psi = phi2 - np.pi/2

        c1 = np.cos(phi)
        c2 = np.cos(theta)
        c3 = np.cos(psi)

        s1 = np.sin(phi)
        s2 = np.sin(theta)
        s3 = np.sin(psi)
        
        R1 = np.array([[c1, -s1, 0],[s1, c1, 0], [0, 0, 1]])
        R2 = np.array([[1, 0, 0],[0, c2, -s2], [0, s2, c2]])
        R3 = np.array([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])

        rot1 = np.dot(R1,position)
        rot2 = np.dot(R2,rot1)
        rot3 = np.dot(R3,rot2)

        return rot3

    def calculate_partial_force(self, position, dx=1e-5):
        # dx = 1e18 # is this causing problems?
        # dx = 1e-5
        fx = partial_derivative(self.potential_function,0,position,dx)
        fy = partial_derivative(self.potential_function,1,position,dx)
        fz = partial_derivative(self.potential_function,2,position,dx)
        
        vec = np.array([-fx,-fy,-fz])
        # ep = 1e-16
        # vec[np.where(np.abs(vec) <= 1000*ep)] = 0
        
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
