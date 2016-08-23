from __future__ import division
import matplotlib.pyplot as plt
import numpy as numpy
import sys
import astropy.units as u
from gala.units import UnitSystem
import gala.potential as gp

from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.potentials.test_spherical_potentials import *
from sidm_orbit_calculation.src.calculation.mass import *
from sidm_orbit_calculation.src.utils.constants import *
import sidm_orbit_calculation.src.potentials.triaxial_BT as BT


potential = 'triaxial_NFW_BT'
host_idx = int(sys.argv[1])
host = HostHalo(idx=host_idx, potential=potential, subs=False)
host.update(host.cosmo.age(0))
# host.q = 0.999
# host.s = 0.99

gravity = GetGravitationalForce(host)

usys = UnitSystem(u.Mpc, u.Gyr, u.Msun, u.degree, u.Mpc/u.Gyr)

M = host.mass_function(host=host, a=0, b=host.R_s)
# v = np.sqrt(G*M/host.R_s)
fg, _ = gravity.calculate_gravitational_force([host.R_s,0,0])
v_c = np.sqrt(fg*host.R_s)

'''
eb = 0.3
ec = 0.5
a = 1
b = np.sqrt(1-eb**2)
c = np.sqrt(1-ec**2)
print a, b, c
'''
a = 1
b = host.q
c = host.s

# print M, v_c, np.sqrt(G*M/host.R_s)

# if host.potential == 'spherical_NFW':
# 	pot = gp.SphericalNFWPotential(v_c = v*u.Mpc/u.Gyr, r_s=R*u.Mpc, units=usys)
# if host.potential == 'triaxial_NFW':
triaxial = gp.LeeSutoTriaxialNFWPotential(v_c =v_c*u.Mpc/u.Gyr, r_s = host.R_s*u.Mpc, a = a, b = b, c = c, units=usys)

x_vals = np.linspace(1e-2*host.R,host.R,100)
y = 0
z = 0

tri = []
tri_gala = []
tri_BT = []
tri_dens = []
tri_dens_gala = []

# ratio = triaxial.density([x_vals[0],y,z]).T[0]/(host.density_function(x_vals[0]*m_to_Mpc,y*m_to_Mpc,z*m_to_Mpc,host)*(m_to_Mpc**3*M_sol))
# print ratio
# ratio = 0.551569218816
# host.rho_s*=ratio
# print host.rho_s

for x in x_vals:
	tri_gala.append(triaxial.value([x,y,z]))
	
	pot = triaxial_NFW_potential(x,y,z,host)
	tri.append(pot)

	pot_BT = BT.triaxial_NFW_potential(x,y,z,host)
	tri_BT.append(pot_BT)

	dens = host.density_function(x,y,z,host)
	tri_dens.append(dens)
	tri_dens_gala.append(triaxial.density([x,y,z]))

tri = np.asarray(tri)
tri_gala = np.asarray(tri_gala)
tri_BT = np.asarray(tri_BT) # /1.33975966654 # density ratio

tri_dens = np.asarray(tri_dens)
tri_dens_gala = np.asarray(tri_dens_gala)

tri_ratio = tri/tri_gala.T[0]
tri_ratio_BT = tri_BT/tri_gala.T[0]
tri_dens_ratio = tri_dens/tri_dens_gala.T[0]

# print tri_dens_ratio.min(), tri_dens_ratio.max()
print tri_ratio_BT.min(), tri_ratio_BT.max()
print tri_ratio.min(), tri_ratio.max()

plt.figure()
plt.plot(x_vals, tri, 'c', linewidth=3, label='triaxial')
plt.plot(x_vals, tri_gala, 'g', linewidth=3, label='gala triaxial')
plt.plot(x_vals, tri_BT, 'b', linewidth=3, label='BT triaxial')
plt.title('potential')
plt.legend()
plt.xscale('log')
plt.grid()

plt.figure()
plt.plot(x_vals, tri_ratio, 'c', linewidth=3, label='triaxial ratio')
plt.plot(x_vals, tri_ratio_BT, 'b', linewidth=3, label='BT triaxial ratio')
plt.xscale('log')
plt.title('potential ratio')
plt.legend()
plt.grid()

'''
plt.figure()
plt.plot(x_vals, tri_dens, 'c--', linewidth=2, label='triaxial')
plt.plot(x_vals, tri_dens_gala, 'g--', linewidth=2, label='gala triaxial')
plt.title('density')
plt.legend()
plt.xscale('log')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.ticklabel_format(style = 'sci', useOffset=False)
ax.plot(x_vals, tri_dens_ratio, 'k', linewidth=2, label='triaxial density ratio')
plt.xscale('log')
plt.title('density ratio')
plt.legend()
'''

# plt.show()


######################################

'''

host_sph = HostHalo(idx=host_idx,potential='spherical_NFW',subs=False)
host_sph.update(host_sph.cosmo.age(0))

M_sph = calculate_spherical_mass(host=host_sph, a=0, b=host_sph.R_s)
v_sph = np.sqrt(G*M_sph/host_sph.R_s)
R_sph = host_sph.R_s

spherical = gp.SphericalNFWPotential(v_c = v_sph*u.Mpc/u.Gyr, r_s=R_sph*u.Mpc, units=usys)

# x = np.linspace(1e20,1e22,100)
# x = np.ones(100)*1e22
# y = np.zeros(100)
# z = np.zeros(100)

# sph_gala = spherical.value(np.vstack([x,y,z]))
# sph = spherical_NFW_potential(x*m_to_Mpc,y*m_to_Mpc,z*m_to_Mpc,host_sph)

sph = []
sph_gala = []
sph_dens = []
sph_dens_gala = []

for x in x_vals:
	sph_gala.append(spherical.value([x,y,z]))
	sph.append(spherical_NFW_potential(x, y, z, host_sph))

	dens = host_sph.density_function(x, y, z, host_sph)
	sph_dens.append(dens)
	sph_dens_gala.append(spherical.density([x,y,z]))

sph = np.asarray(sph)
sph_gala = np.asarray(sph_gala)
sph_dens = np.asarray(sph_dens)
sph_dens_gala = np.asarray(sph_dens_gala)

sph_ratio = sph / sph_gala.T[0]
sph_dens_ratio = sph_dens/sph_dens_gala.T[0]

print sph_dens_ratio.min(), sph_dens_ratio.max()

plt.figure()
plt.plot(x_vals, sph_dens, 'c', lw=2, label='spherical')
plt.plot(x_vals, sph_dens_gala, 'g--', lw=2, label='spherical gala')
plt.title('spherical NFW density profile')
plt.legend()
plt.xscale('log')

plt.figure()
plt.plot(x_vals, sph_dens_ratio, 'k', linewidth=2, label='spherical density ratio')
plt.xscale('log')
plt.title('density ratio')
plt.legend()

plt.figure()
plt.plot(x_vals, sph, 'b', linewidth=2, label='spherical')
plt.plot(x_vals, sph_gala, 'k--', linewidth=2, label='gala spherical')
plt.legend()
plt.xscale('log')

plt.figure()
plt.plot(x_vals, sph_ratio, 'b', linewidth=2, label='spherical ratio')
plt.xscale('log')
plt.legend()
'''

plt.show()
