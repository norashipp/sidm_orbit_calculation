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

'''
tri = []
tri_gala = []

host = HostHalo(M=1e14,potential='triaxial_NFW')

host.s = float(sys.argv[1])
host.q = float(sys.argv[2])

usys = UnitSystem(u.meter, u.second, u.kilogram, u.degree, u.meter/u.second)

M_tri = calculate_triaxial_mass(host=host, a=0, b=host.R_s)
v_tri = np.sqrt(G*M_tri/host.R_s)*(s_to_Gyr/m_to_Mpc)
R = host.R_s/m_to_Mpc

triaxial = gp.LeeSutoTriaxialNFWPotential(v_c =v_tri*u.m/u.s, r_s = R*u.m, a = 1, b = host.q, c = host.s, units=usys)

x_vals = np.linspace(1e20,1e22,100)
y = 0
z = 0

for x in x_vals:
	tri_gala.append(triaxial.value([x,y,z]))

x_vals *= m_to_Mpc
for x in x_vals:
	pot = triaxial_NFW_potential(x,y,z,host)
	tri.append(pot)

x_vals /= m_to_Mpc

tri = np.asarray(tri)
tri_gala = np.asarray(tri_gala)
tri_ratio = tri*(s_to_Gyr/m_to_Mpc)**2/tri_gala.T[0]
tri *= (s_to_Gyr/m_to_Mpc)**2

plt.figure()
plt.plot(x_vals, tri, 'c--', linewidth=2, label='triaxial')
plt.plot(x_vals, tri_gala, 'g--', linewidth=2, label='gala triaxial')
plt.legend()
plt.xscale('log')

plt.figure()
plt.plot(x_vals, tri_ratio, 'k', linewidth=2, label='triaxial ratio')
plt.xscale('log')
plt.legend()

# plt.show()
'''

######################################

sph = []
sph_gala = []

host_sph = HostHalo(M=1e14,potential='spherical_NFW')

usys = UnitSystem(u.meter, u.second, u.kilogram, u.degree, u.meter/u.second)

M_sph = calculate_spherical_mass(host=host_sph, a=0, b=host_sph.R_s)
v_sph = np.sqrt(G*M_sph/host_sph.R_s)*(s_to_Gyr/m_to_Mpc)
R_sph = host_sph.R_s/m_to_Mpc

spherical = gp.SphericalNFWPotential(v_c = v_sph*u.m/u.s, r_s=R_sph*u.m, units=usys)

x = np.linspace(1e20,1e22,100)
# x = np.ones(100)*1e22
y = np.zeros(100)
z = np.zeros(100)

# sph_gala = spherical.value(np.vstack([x,y,z]))
# sph = spherical_NFW_potential(x*m_to_Mpc,y*m_to_Mpc,z*m_to_Mpc,host_sph)

for i in range(len(x)):
	sph_gala.append(spherical.value([x[i],y[i],z[i]]))
	sph.append(spherical_NFW_potential(x[i]*m_to_Mpc,y[i]*m_to_Mpc,z[i]*m_to_Mpc,host_sph))

sph = np.asarray(sph)
sph_gala = np.asarray(sph_gala)

sph_gala *= (m_to_Mpc / s_to_Gyr) ** 2
sph_ratio = sph / sph_gala.T[0]

plt.figure()
plt.plot(x, sph, 'b', linewidth=2, label='spherical')
plt.plot(x, sph_gala, 'k--', linewidth=2, label='gala spherical')
plt.legend()
plt.xscale('log')

plt.figure()
plt.plot(x, sph_ratio, 'b', linewidth=2, label='spherical ratio')
plt.xscale('log')
plt.legend()

plt.show()
