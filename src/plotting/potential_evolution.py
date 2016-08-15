import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from gala.units import UnitSystem
import gala.potential as gp

from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.utils.constants import *

dpi = 175
fontsize = 12
plt.rc('savefig', dpi=dpi)
plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.rc('xtick.major', pad=5)
plt.rc('xtick.minor', pad=5)
plt.rc('ytick.major', pad=5)
plt.rc('ytick.minor', pad=5)

host_halo_mass = 1e14
potential = 'spherical_NFW'
host_idx = 40
tmax = 10
dt = 1

host = HostHalo(M=host_halo_mass, potential=potential, idx=host_idx, tmax=tmax)

usys = UnitSystem(u.meter, u.second, u.kilogram, u.degree, u.meter/u.second)

def gala_spherical(pos):
	M = host.mass_function(host=host, a=0, b=host.R_s)
	v = np.sqrt(G*M/host.R_s)*(s_to_Gyr/m_to_Mpc)
	R = host.R_s/m_to_Mpc
	spherical = gp.SphericalNFWPotential(v_c = v*u.m/u.s, r_s=R*u.m, units=usys)
	return spherical.value(pos)

time = host.cosmo.age(0) - tmax

xx = np.linspace(0.01,1,100)*host.R
yy = np.zeros(100)
zz = np.zeros(100)

pos = np.vstack([xx,yy,zz])/m_to_Mpc

plt.figure()

while time < tmax:
	time+=dt

	z = host.cosmo.age(time,inverse=True)
	print 'z = ', z

	host.update(time)

	pp = host.potential_function(xx,yy,zz,host)
	pg = gala_spherical(pos)*(m_to_Mpc/s_to_Gyr)**2

	plt.plot(xx/host.R, pp, label='z = %.2f' %z, lw=2)
	plt.plot(xx/host.R, pg, '--', label='z = %.2f' %z, lw=2)

plt.legend(loc='lower right')
plt.xlabel(r'$\mathrm{x/R_h}$')
plt.ylabel(r'$\mathrm{potential\ (Mpc^2/Gyr^2)}$')
plt.xscale('log')
plt.show()
