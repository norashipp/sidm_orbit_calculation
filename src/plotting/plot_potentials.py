import matplotlib.pyplot as plt
import numpy as numpy
import sys
import astropy.units as u
from gala.units import UnitSystem
import gala.potential as gp

from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.potentials.test_spherical_potentials import *

tri = []
sph = []
tri_gala = []
sph_gala = []

host = HostHalo(M=1e13,potential='triaxial_NFW')
host_sph = HostHalo(M=1e13,potential='spherical_NFW')

host.s = float(sys.argv[1])
host.q = float(sys.argv[2])

usys = UnitSystem(u.meter, u.second, u.kilogram, u.degree, u.meter/u.second)

triaxial = gp.LeeSutoTriaxialNFWPotential(v_c = host.v*u.m/u.s, r_s = host.R_s*u.m, a = host.R*u.m, b = host.R*host.q*u.m, c = host.R*host.s*u.m, units=usys)
spherical = gp.SphericalNFWPotential(v_c = host_sph.v*u.m/u.s, r_s=host.R_s*u.m, units=usys)

x_vals = np.linspace(1e20,1e22,100)
y = 0
z = 0

for x in x_vals:
	tri.append(triaxial_NFW_potential(x,y,z,host))
	sph.append(spherical_NFW_potential(x,y,z,host_sph))
	tri_gala.append(triaxial.value([x,y,z]))
	sph_gala.append(spherical.value([x,y,z]))

tri = np.asarray(tri)
tri_gala = np.asarray(tri_gala)
tri_ratio = tri/tri_gala.T[0]
# print tri_ratio

plt.figure()
plt.plot(x_vals, tri, 'c--', linewidth=2, label='triaxial')
plt.plot(x_vals, sph, 'b', linewidth=2, label='spherical')
plt.plot(x_vals, tri_gala, 'g--', linewidth=2, label='gala triaxial')
plt.plot(x_vals,sph_gala, 'k', linewidth=2, label='gala spherical')
plt.legend()
plt.xscale('log')

plt.figure()
plt.plot(x_vals, tri_ratio, 'k', linewidth=2, label='triaxial ratio')
plt.xscale('log')
plt.legend()

plt.show()