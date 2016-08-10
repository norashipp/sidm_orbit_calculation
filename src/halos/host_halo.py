import numpy as np
from scipy.integrate import quad

from colossus.cosmology import cosmology
from colossus.halo.mass_so import M_to_R
from colossus.halo.concentration import concentration
from colossus.halo.mass_defs import changeMassDefinition

from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.utils.constants import *
import sidm_orbit_calculation.src.potentials.density as density
import sidm_orbit_calculation.src.calculation.mass as mass

class HostHalo:

    def __init__(self,M,potential,s=0.39,q=0.58):
        my_cosmo = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.045714, 'sigma8': 0.82, 'ns': 0.96}
        self.cosmo = cosmology.setCosmology('my_cosmo', my_cosmo)

        self.potential = potential

        self.density_function = density.density_dict[self.potential]
        self.mass_function = mass.mass_dict[self.potential]

        self.M = M # work in units of solar mass, virial mass
        
        self.z = 0.0 # from merger trees?
        
        self.s = s # axis ratios
        self.q = q

        self.c = self.concentration()
        self.R = self.virial_radius()
        self.v = self.virial_velocity()

        self.R_s = self.scale_radius()
        
        self.rho_s = 1
        self.rho_s = self.scale_density()
        
        self.rho = 0

        self.density_array = []

    def virial_mass(self):
        M, R, c = changeMassDefinition(M_200m, c_200m, z, '200m', 'vir')

    def concentration(self):
        c = concentration(self.M, 'vir', self.z, model='diemer15')
        return c

    def virial_radius(self):
        radius = M_to_R(self.M,self.z,'vir') # correct mdef? # returns kpc
        return radius/1000 # Mpc # /m_to_kpc

    def scale_radius(self):
        return self.R/self.c # R_vir / c_vir

    def virial_velocity(self):
        # velocity = np.sqrt(2*G*self.M*M_sol/self.R)
        velocity = np.sqrt(2*G*self.M/self.R)
        return velocity

    def scale_density(self):
        return self.M/self.mass_function(self,0,self.R)
