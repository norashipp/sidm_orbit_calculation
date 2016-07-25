import numpy as np

from colossus.cosmology import cosmology
from colossus.halo.mass_so import M_to_R
from colossus.halo.concentration import concentration

from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.utils.constants import *

class HostHalo:

    def __init__(self,M,potential,s=0.39,q=0.58):
        my_cosmo = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.045714, 'sigma8': 0.82, 'ns': 0.96}
        self.cosmo = cosmology.setCosmology('my_cosmo', my_cosmo)

        self.potential = potential
        
        self.M = M # work in units of solar mass, virial mass
        
        self.z = 0.0 # from merger trees?
        
        self.s = s # axis ratios
        self.q = q

        self.c = self.concentration()
        self.R = self.virial_radius()
        self.v = self.virial_velocity()

        self.M *= M_sol # need to go through and edit all files to change units

        self.R_s = self.scale_radius()
        self.rho_s = self.scale_density()

        self.update_density([self.R,0])
        self.density_array = [self.rho]

    def update_density(self,position):
        # for now just NFW - should have a dictionary like for potentials
        r = np.sqrt(position[0]**2 + position[1]**2)
        x = r/self.R_s
        self.rho = self.rho_s/(x*(1+x)**2)
        # self.rho = 1e-40 # average - for analytical comparison
        
    def concentration(self):
        c = concentration(self.M, 'vir', self.z, model='diemer15')
        return c

    def virial_radius(self):
        radius = M_to_R(self.M,self.z,'vir') # correct mdef? # returns kpc
        return radius/m_to_kpc

    def scale_radius(self):
        return self.R/self.c # R_vir / c_vir correct

    def virial_velocity(self):
        velocity = np.sqrt(2*G*self.M*M_sol/self.R)
        return velocity

    def scale_density(self):
        # kg/m3
        return self.M/(4*np.pi*self.R_s**3*(np.log(1+self.c)-self.c/(1+self.c)))
        