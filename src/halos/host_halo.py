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

    def __init__(self, M, potential, s = 0.99, q = 0.999, a = 1, R_s = None):
        # input from merger tree: R_s, a, M_200m, b_to_a, c_to_a

        my_cosmo = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.045714, 'sigma8': 0.82, 'ns': 0.96}
        self.cosmo = cosmology.setCosmology('my_cosmo', my_cosmo)

        self.potential = potential

        self.density_function = density.density_dict[self.potential]
        self.mass_function = mass.mass_dict[self.potential]

        # axis ratios, q > s
        self.s = s # c_to_a
        self.q = q # b_to_a

        self.z = self.redshift(a)

        self.R = self.virial_radius(M)

        self.c = self.concentration(M, '200m')

        self.M = M # virial_mass(M) # work in 200m?
        
        self.v = self.virial_velocity()

        self.R_s = self.scale_radius() # this is in the merger tree's...check consistency
        print 'scale radius'
        print self.R_s
        print R_s
        print 
        
        self.rho_s = 1
        self.rho_s = self.scale_density()
        
        self.rho = 0

        self.density_array = []

    def redshift(self, a):
        return 1/a - 1

    def virial_radius(self, M):
        radius = M_to_R(M, self.z, '200m') # correct mdef? # returns kpc
        return radius/1000 # Mpc

    def concentration(self, M, mdef='200m'):
            c = concentration(M, mdef, self.z, model='diemer15')
            return c

    # def virial_mass(self, M_200m):
    #     # unnecessary?
    #     M, R, c = changeMassDefinition(M_200m, c_200m, z, '200m', 'vir')
    #     return M

    def scale_radius(self):
        return self.R/self.c # R_vir / c_vir

    def virial_velocity(self):
        return np.sqrt(2*G*self.M/self.R)
        
    def scale_density(self):
        return self.M/self.mass_function(self,0,self.R)
