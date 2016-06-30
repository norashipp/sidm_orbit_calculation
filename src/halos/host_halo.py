import numpy as np

from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.utils.constants import *

class HostHalo:

    def __init__(self,M,potential):
        self.M = M
        self.potential = potential
        self.c_vir = self.concentration()
        self.R_200 = self.virial_radius()
        self.R_s = self.scale_radius()
        self.v_200 = self.virial_velocity()
        self.rho_s = self.scale_density()

        self.update_density([self.R_200,0])
        self.density_array = [self.rho]

    def update_density(self,position):
        # for now just NFW
        r = np.sqrt(position[0]**2 + position[1]**2)
        x = r/self.R_s
        self.rho = self.rho_s/(x*(1+x)**2)
        self.rho = 1e-40 # average
        
    def concentration(self):
        c_vir = 6.0 # typical for M~1e14, from Benedikt's paper, fig 2?
        return c_vir

    def virial_radius(self):
        radius = 1.e22 # meters, ~10 Mpc
        return radius

    def scale_radius(self):
        return self.R_200/self.c_vir

    def virial_velocity(self):
        velocity = 1.e6 # m/s
        return velocity

    def scale_density(self):
        # kg/m3
        return self.M/(4*np.pi*self.R_s**3*(np.log(1+self.c_vir)-self.c_vir/(1+self.c_vir)))
        