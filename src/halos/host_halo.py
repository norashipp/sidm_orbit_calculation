import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

from colossus.cosmology import cosmology
from colossus.halo.mass_so import M_to_R
from colossus.halo.concentration import concentration
from colossus.halo.mass_defs import changeMassDefinition

from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.utils.constants import *
import sidm_orbit_calculation.src.potentials.density as density
import sidm_orbit_calculation.src.calculation.mass as mass
from sidm_orbit_calculation.src.utils.setup import *
import sidm_orbit_calculation.src.merger_tree.cluster as cluster
import sidm_orbit_calculation.src.potentials.test_spherical_potentials as potentials

HOMEDIR = home_directory()

class HostHalo:

    def __init__(self, M, potential, idx=None, tmax=None, R_s=None, s=0.99, q=0.999, a=1):
        # input from merger tree: R_s, a, M_200m, b_to_a, c_to_a

        my_cosmo = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.045714, 'sigma8': 0.82, 'ns': 0.96}
        self.cosmo = cosmology.setCosmology('my_cosmo', my_cosmo)

        self.host_idx = idx

        self.potential = potential

        self.potential_function = potentials.potential_dict[potential]
        self.density_function = density.density_dict[self.potential]
        self.mass_function = mass.mass_dict[self.potential]

        if idx != None:
            self.merger_tree(idx)
            t0 = self.cosmo.age(0)-tmax
            print self.cosmo.age(0)
            print tmax
            print 'initial time = ', t0
            self.update(t0)

            self.subs = SubHalos(HOMEDIR + 'subs/sub_%d.dat' % self.host_idx)

        else:
            # axis ratios, q > s
            self.s = s # c_to_a
            self.q = q # b_to_a

            self.z = self.redshift(a)

            self.M = M # virial_mass(M) # work in 200m?

            self.R = self.virial_radius()

            self.R_s = R_s
            
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

    def virial_radius(self):
        radius = M_to_R(self.M, self.z, '200m') # correct mdef? # returns kpc
        return radius/1000 # Mpc

    def concentration(self):
        if self.R_s:
            return self.R/self.R_s
        else: 
            c = concentration(self.M, '200m', self.z, model='diemer15')
            self.R_s = self.R/c
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

    def merger_tree(self,host_idx=40):
        print 'Importing host halo parameters from merger tree...'
        # homedir = home_directory()
        self.hs = cluster.HostHalos(HOMEDIR + 'sidm_orbit_calculation/src/merger_tree/clusters.dat')
        
        aa = self.hs[self.host_idx].a 
        zz = self.redshift(aa)
        tt = self.cosmo.age(zz) # Gyr
        print 'time range = ', tt.min(), tt.max()

        self.M_sp = interp1d(tt,self.hs[self.host_idx].m_200m)
        self.R_s_sp = interp1d(tt,self.hs[self.host_idx].r_s)
        self.q_sp = interp1d(tt,self.hs[self.host_idx].b_to_a)
        self.s_sp = interp1d(tt,self.hs[self.host_idx].c_to_a)

    def update(self, time):
        self.M = self.M_sp(time)
        self.z = self.cosmo.age(time,inverse=True)
        self.R_s = self.R_s_sp(time)
        self.q = self.q_sp(time)
        self.s = self.s_sp(time)

        self.R = self.virial_radius()
        self.c = self.concentration()        
        self.v = self.virial_velocity()
