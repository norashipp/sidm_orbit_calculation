from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from sidm_orbit_calculation.src.merger_tree.cluster import *
from sidm_orbit_calculation.src.utils.setup import *
from sidm_orbit_calculation.src.utils.geometry import *
hs = HostHalos(HOMEDIR + 'sidm_orbit_calculation/src/merger_tree/clusters.dat')

# 51 hosts

a_theta = []
a_phi = []
mass = []
scale_radius = []
b_to_a = []
c_to_a = []

for i in range(len(hs)):
	host = hs[i]
	axis = np.array([host.ax[-1],host.ay[-1],host.az[-1]])
	_, theta, phi = spherical_coordinates(axis)
	a_theta.append(theta)
	a_phi.append(phi)
	mass.append(host.m_200m[-1])
	scale_radius.append(host.r_s[-1])
	b_to_a.append(host.b_to_a[-1])
	c_to_a.append(host.c_to_a[-1])

nbins = 10

plt.figure()
plt.hist(np.asarray(a_theta)/np.pi,bins=nbins,histtype='step',color='c',normed=False,lw=3,label='theta')
plt.hist(np.asarray(a_phi)/np.pi,bins=nbins,histtype='step',color='g',normed=False,lw=3,label='phi')
plt.legend()
plt.xlabel('Pi Radians')
plt.title('Direction of Major Axis')
plt.savefig('host_halo_major_axes.png')

plt.figure()
plt.hist(mass,bins=nbins,histtype='step',color='c',normed=False,lw=3)
plt.xlabel('M_200m (M_sun)')
plt.title('Mass')
plt.savefig('host_halo_masses.png')

plt.figure()
plt.hist(scale_radius,bins=nbins,histtype='step',color='c',normed=False,lw=3)
plt.xlabel('r_s (Mpc)')
plt.title('Scale Radius')
plt.savefig('host_halo_scale_radii.png')

plt.figure()
plt.hist(b_to_a,bins=nbins,histtype='step',color='c',normed=False,lw=3,label='b_to_a')
plt.hist(c_to_a,bins=nbins,histtype='step',color='g',normed=False,lw=3,label='c_to_a')
plt.legend()
plt.title('Axis Ratios')
plt.savefig('host_halo_axis_ratios.png')

# plt.show()