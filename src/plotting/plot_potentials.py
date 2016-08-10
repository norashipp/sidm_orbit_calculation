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

tri = []
sph = []
tri_gala = []
sph_gala = []

host = HostHalo(M=1e14,potential='triaxial_NFW')
host_sph = HostHalo(M=1e14,potential='spherical_NFW')

host.s = float(sys.argv[1])
host.q = float(sys.argv[2])

usys = UnitSystem(u.meter, u.second, u.kilogram, u.degree, u.meter/u.second)

M_tri = calculate_triaxial_mass(host=host, a=0, b=host.R_s)
# print '1 ', M_tri
# print '2 ', host.M
M_sph = calculate_spherical_mass(host=host_sph, a=0, b=host_sph.R_s)
# M_tri = 4.54e+11*M_sol

v_tri = np.sqrt(G*M_tri/host.R_s)*(s_to_Gyr/m_to_Mpc)
v_sph = np.sqrt(G*M_sph/host_sph.R_s)*(s_to_Gyr/m_to_Mpc)

R = host.R_s/m_to_Mpc
R_sph = host_sph.R_s/m_to_Mpc

triaxial = gp.LeeSutoTriaxialNFWPotential(v_c =v_tri*u.m/u.s, r_s = R*u.m, a = 1, b = host.q, c = host.s, units=usys)
spherical = gp.SphericalNFWPotential(v_c = v_sph*u.m/u.s, r_s=R_sph*u.m, units=usys)

x_vals = np.linspace(1e20,1e22,100)
y = 0
z = 0

for x in x_vals:
	# tri.append(triaxial_NFW_potential(x,y,z,host))
	# sph.append(spherical_NFW_potential(x,y,z,host_sph))
	tri_gala.append(triaxial.value([x,y,z]))
	sph_gala.append(spherical.value([x,y,z]))

x_vals *= m_to_Mpc
for x in x_vals:
	pot = triaxial_NFW_potential(x,y,z,host)
	# print x, pot
	# pot *= (m_to_Mpc/s_to_Gyr) ** -2
	tri.append(pot)
	sph.append(spherical_NFW_potential(x,y,z,host_sph))

x_vals /= m_to_Mpc

tri = np.asarray(tri)
tri_gala = np.asarray(tri_gala)
sph = np.asarray(sph)
sph_gala = np.asarray(sph_gala)
tri_ratio = tri*(s_to_Gyr/m_to_Mpc)**2/tri_gala.T[0]
sph_ratio = sph*(s_to_Gyr/m_to_Mpc)**2/sph_gala.T[0]

tri *= (s_to_Gyr/m_to_Mpc)**2

plt.figure()
plt.plot(x_vals, tri, 'c--', linewidth=2, label='triaxial')
# plt.plot(x_vals, sph, 'b', linewidth=2, label='spherical')
plt.plot(x_vals, tri_gala, 'g--', linewidth=2, label='gala triaxial')
# plt.plot(x_vals,sph_gala, 'k', linewidth=2, label='gala spherical')
plt.legend()
plt.xscale('log')

plt.figure()
plt.plot(x_vals, tri_ratio, 'k', linewidth=2, label='triaxial ratio')
# plt.plot(x_vals, sph_ratio, 'b', linewidth=2, label='spherical ratio')
plt.xscale('log')
plt.legend()

plt.show()
