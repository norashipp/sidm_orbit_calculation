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

    def __init__(self, M, potential, R_s = None, s = 0.99, q = 0.999, a = 1):
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

        self.R_s = R_s

        self.M = M # virial_mass(M) # work in 200m?
        
        self.c = self.concentration()
        
        self.v = self.virial_velocity()

        self.rho_s = 1
        self.rho_s = self.scale_density()
        
        self.rho = 0

        self.density_array = []

        print 
        print 'HOST HALO'
        print 'axis ratios = %.3f, %.3f' %(self.s, self.q)
        print 'redshift = %.3f' %self.z
        print 'mass = %.3e' %self.M
        print 'virial radius = %.3e' %self.R
        print 'scale radius = %.3e' %self.R_s
        print 'concentration = %.3f' %self.c
        print 'virial velocity = %.3e' %self.v
        print 'scale density = %.3e' %self.rho_s
        print 

    def redshift(self, a):
        return 1/a - 1

    def virial_radius(self, M):
        radius = M_to_R(M, self.z, '200m') # correct mdef? # returns kpc
        return radius/1000 # Mpc

    def concentration(self):
        if self.R_s:
            return self.R/self.R_s
        else: 
            c = concentration(self.M, '200m', self.z, model='diemer15')
            self.R_s = self.R/self.c
            return c

    def virial_velocity(self):
        return np.sqrt(2*G*self.M/self.R)
        
    def scale_density(self):
        return self.M/self.mass_function(self,0,self.R)

    # def virial_mass(self, M_200m):
    #     # unnecessary?
    #     M, R, c = changeMassDefinition(M_200m, c_200m, z, '200m', 'vir')
    #     return M

    # def scale_radius(self):
    #     return self.R/self.c # R_vir / c_vir
    