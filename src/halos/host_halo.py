from __future__ import division
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from time import sleep
from scipy.interpolate import UnivariateSpline

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
from sidm_orbit_calculation.src.halos.subhalo import *
from sidm_orbit_calculation.src.utils.geometry import *

class HostHalo:

    # def __init__(self, M, potential, idx=None, tmax=None, R_s=None, s=0.99, q=0.999, a=1):
    def __init__(self, idx, potential, sigma=0, subs=True, scale_density=True):
        # input from merger tree: R_s, a, M_200m, b_to_a, c_to_a

        my_cosmo = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.0469, 'sigma8': 0.82, 'ns': 0.95}
        self.cosmo = cosmology.setCosmology('my_cosmo', my_cosmo)

        self.host_idx = idx
        self.sigma = sigma
        self.potential = potential

        self.potential_function = potentials.potential_dict[potential]
        self.density_function = density.density_dict[self.potential]
        self.mass_function = mass.mass_dict[self.potential]
        
        self.initiate_host(scale_density=scale_density)
        if subs: self.initiate_subhalos()
        
        self.rho = 0
        self.density_array = []

        '''
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
        '''

        '''
        if idx != None:
            self.initiate_host()
            t0 = self.cosmo.age(0)-tmax
            # print self.cosmo.age(0)
            # print tmax
            # print 'initial time = ', t0
            self.update(t0)

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
        '''

    def redshift(self, a):
        return 1/a - 1

    def virial_radius(self, M, z):
        radius = M_to_R(M*self.cosmo.h, z, '200m')/self.cosmo.h # returns kpc
        return radius/1000 # Mpc

    def concentration(self):
        return self.R/self.R_s
        # if self.R_s:
        #     return self.R/self.R_s
        # else: 
        #     c = concentration(self.M, '200m', self.z, model='diemer15')
        #     self.R_s = self.R/c
        #     return c

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

    def initiate_host(self,scale_density):
        print 'Importing host %i parameters from merger tree...' %self.host_idx
        hs = cluster.HostHalos(HOMEDIR + 'sidm_orbit_calculation/src/merger_tree/clusters.dat')
        aa = hs[self.host_idx].a 
        zz = self.redshift(aa)
        tt = self.cosmo.age(zz) # Gyr
        # print 'time range = ', tt.min(), tt.max()

        self.M_sp = interp1d(tt,hs[self.host_idx].m_200m/self.cosmo.h)
        self.R_s_sp = interp1d(tt,hs[self.host_idx].r_s/(self.cosmo.h*(1+zz)))
        self.q_sp = interp1d(tt,hs[self.host_idx].b_to_a)
        self.s_sp = interp1d(tt,hs[self.host_idx].c_to_a)
        
        self.ax_sp = interp1d(tt,hs[self.host_idx].ax)
        self.ay_sp = interp1d(tt,hs[self.host_idx].ay)
        self.az_sp = interp1d(tt,hs[self.host_idx].az)

        if scale_density:
            rho = []
            self.rho_s = 1
            for t in tt:
                self.z = self.cosmo.age(t,inverse=True)
                self.M = self.M_sp(t)
                self.R_s = self.R_s_sp(t)
                self.R = self.virial_radius(self.M, self.z)
                self.q = self.q_sp(t)
                self.s = self.s_sp(t)
                rho.append(self.scale_density())

            rho = np.asarray(rho)
            lrho = np.log(rho)
        else:
            print 'skipping scale density calculation'
            lrho = np.ones_like(tt)

        self.lrho_s_sp = interp1d(tt,lrho)


    def initiate_subhalos(self):
        subs = cluster.SubHalos(HOMEDIR + 'sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat' % self.host_idx)
        print "Initiating subhahlos..."
        self.subhalos = []
        skip = 0
        sub_count = 0
        # for sub in subs:
        for i in range(len(subs)-1):
            sub = subs[i]
            sub_count+=1

            vthresh = 70 # km/s
            if sub.v_max[-1] < vthresh:
                self.subhalos.append(None)
                skip+=1
                continue

            aa = sub.a
            zz = self.redshift(aa)
            tt = self.cosmo.age(zz) # Gyr
            
            if len(tt) < 10:
                self.subhalos.append(None)
                continue

            hostM = self.M_sp(tt)
            hostR = self.virial_radius(hostM, zz)
            dd = np.sqrt(sub.rel_x*sub.rel_x + sub.rel_y*sub.rel_y + sub.rel_z*sub.rel_z) # Mpc
            ratio_to_tt = interp1d(dd/hostR,tt)            
            
            # t0 = tt[0]
            # t0 = tt[-10]
            try:
                t0 = ratio_to_tt(1) # determine when subhalo is at a distance of n * host.R
            except:
                self.subhalos.append(None)
                skip+=1
                continue
            
            x = sub.rel_x/(self.cosmo.h*(1/sub.a))
            y = sub.rel_y/(self.cosmo.h*(1/sub.a))
            z = sub.rel_z/(self.cosmo.h*(1/sub.a))

            t_to_x = interp1d(tt,x)
            t_to_y = interp1d(tt,y)
            t_to_z = interp1d(tt,z)
            # t_to_vx = interp1d(tt,sub.rel_vx*1000*m_to_Mpc/s_to_Gyr)
            # t_to_vy = interp1d(tt,sub.rel_vy*1000*m_to_Mpc/s_to_Gyr)
            # t_to_vz = interp1d(tt,sub.rel_vz*1000*m_to_Mpc/s_to_Gyr)
            t_to_m = interp1d(tt,sub.m_200m/self.cosmo.h)

            x0, y0, z0 = t_to_x(t0), t_to_y(t0), t_to_z(t0)
            # vx0, vy0, vz0 = t_to_vx(t0), t_to_vy(t0), t_to_vz(t0)
            m0 = t_to_m(t0)

            # print 'checking lengths: ', len(tt), len(sub.rel_x/(self.cosmo.h*(1/sub.a)))
            s = 0.05
            xsp = UnivariateSpline(tt,x,s=s)
            ysp = UnivariateSpline(tt,y,s=s)
            zsp = UnivariateSpline(tt,z,s=s)

            vxsp = xsp.derivative()
            vysp = ysp.derivative()
            vzsp = zsp.derivative()

            '''
            plt.figure(figsize=(5,5))
            plt.plot(tt,xsp(tt))
            plt.plot(tt,ysp(tt))
            plt.plot(tt,zsp(tt))

            plt.figure(figsize=(5,5))
            plt.plot(tt,vxsp(tt))
            plt.plot(tt,vysp(tt))
            plt.plot(tt,vzsp(tt))

            plt.show()
            '''

            vx0 = vxsp(t0)
            vy0 = vysp(t0)
            vz0 = vzsp(t0)

            initial_position = np.array([x0, y0, z0])
            initial_momentum = np.array([vx0, vy0, vz0])

            # print 'initial parameters'
            # print initial_position
            # print initial_momentum

            # initial_position = rotate(initial_position)
            # initial_momentum = rotate(initial_momentum)

            subhalo = Subhalo(host=self, M=m0, initial_position=initial_position, initial_momentum=initial_momentum, t0=t0, mass_spline=t_to_m)
            self.subhalos.append(subhalo)

            '''
            # print i
            if i == 33:
                print t0, initial_position, initial_momentum
                print t0, t_to_x(t0)
                times = np.linspace(tt.min(),tt.max(),1000)
                plt.figure()
                plt.plot(tt,sub.rel_x,'b.',markersize=5)
                plt.plot(times,t_to_x(times),'c',lw=2)
                plt.plot(t0,t_to_x(t0),'c*')
                plt.show()
            '''

            '''
            print 'subhalo'
            print 't0 = ', t0
            print 'x0 = ', x0
            print 'p0 = ', p0
            '''

            '''
            print 'kept: ', len(self.subhalos)
            print 'skipped: ', skip
            print 'total: ', sub_count
            print 
            '''

    def rotate_orbit(velocity):
        u = uniform(0,1)
        v = uniform(0,1)

        phi = 2*np.pi*u
        theta = np.arccos(2*v - 1)
        
        # cos_theta = uniform(-1,1)
        # phi = uniform(0,2*np.pi)

        # sin_theta = np.sqrt(1-cos_theta**2)
        north_south = choice([-1,1])
        if north_south == -1: theta+=np.pi
        # sin_theta *= north_south

        # rotate around x by theta
        # rotate around z by phi
        # cos theta uniformly distributed
        # phi uniformly distributed
        # randomly select between northern and southern hemisphere

        Rx = np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
        Rz = np.array([[np.cos(phi), -np.sin(phi), 0], [np.sin(phi), np.cos(phi), 0], [0, 0, 1]])

        velocity_rotated = np.dot(Rx,velocity)
        velocity_rotated = np.dot(Rz, velocity_rotated)

        position = np.array([1, 0, 0])
        position_rotated = np.dot(Rx, position)
        position_rotated = np.dot(Rz, position_rotated)

        return position_rotated, velocity_rotated

    def update(self, time):
        self.M = self.M_sp(time)
        self.z = self.cosmo.age(time,inverse=True)
        self.R_s = self.R_s_sp(time)
        self.q = self.q_sp(time)
        self.s = self.s_sp(time)

        self.R = self.virial_radius(self.M, self.z)
        self.c = self.concentration()        
        self.v = self.virial_velocity()

        self.rho_s = np.exp(self.lrho_s_sp(time))

        time = self.cosmo.age(0) # maybe evolving axis direction is causing weirdness
        self.ax = self.ax_sp(time)
        self.ay = self.ay_sp(time)
        self.az = self.az_sp(time)

        '''
        print 
        print 'HOST HALO'
        print 'redshift = %.3f' %self.z
        print 'mass = %.3e' %self.M
        print 'axis ratios = %.3f, %.3f' %(self.s, self.q)
        print 'virial radius = %.3e' %self.R
        print 'scale radius = %.3e' %self.R_s
        print 'concentration = %.3f' %self.c
        print 'virial velocity = %.3e' %self.v
        print 'scale density = %.3e' %self.rho_s
        print 
        # sleep(2)
        '''
