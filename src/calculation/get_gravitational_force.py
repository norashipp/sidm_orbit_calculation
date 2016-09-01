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
        axis = np.array([self.host.ax,self.host.ay,self.host.az])
        
        r = linalg.norm(axis)

        costh = self.host.az/r
        sinth = np.sqrt(1 - costh**2)
        
        # rxy = linalg.norm(axis[:2])
        # cosph = self.host.ax/rxy
        # sinph = self.host.ay/rxy
        
        cosph = self.host.ax/(r*sinth)
        sinph = self.host.ay/(r*sinth)

        print 'phi = ', np.arccos(cosph), np.arcsin(sinph)
        print 'theta = ', np.arccos(costh), np.arcsin(sinth)

        # cospsi = 1
        # sinpsi = 0

        # r = np.sqrt(self.host.ax*self.host.ax + self.host.ay*self.host.ay + self.host.az*self.host.az)
        # theta = np.arccos(self.host.az/r)
        # phi = np.arctan2(self.host.ay,self.host.ax)

        
        phi = np.pi/4
        theta = np.pi/2

        # r, theta, phi = spherical_coordinates(axis)
        # print r, theta, phi
        cosph = np.cos(phi)
        sinph = np.sin(phi)
        costh = np.cos(theta)
        sinth = np.sin(theta)
        
        a11 = cosph
        a12 = sinph
        a13 = 0
        a21 = -costh*sinph
        a22 = costh*cosph
        a23 = sinth
        a31 = sinth*sinph
        a32 = -sinth*cosph
        a33 = costh

        A = np.array([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]])

        return np.dot(A,position)

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
