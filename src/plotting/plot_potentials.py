import matplotlib.pyplot as plt
import numpy as numpy
import sys
import astropy.units as u
from gala.units import UnitSystem
import gala.potential as gp

from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.potentials.test_spherical_potentials import *


x_vals = np.linspace(1e18,1e22,100)

tri = []
# sph = []
tri_gala = []

host = HostHalo(M=1e13,potential='triaxial_NFW')
# host_sph = HostHalo(M=1e13,potential='spherical_NFW')

host.s = float(sys.argv[1])
host.q = float(sys.argv[2])

print host.v
print host.s
print host.q
print host.R
print host.R*host.s
print host.R*host.q

usys = UnitSystem(u.meter, u.second, u.kilogram, u.degree, u.meter/u.second)

triaxial = gp.LeeSutoTriaxialNFWPotential(host.v,host.R_s,host.R,host.R*host.q,host.R*host.s,units=usys)


y = 0
z = 0

for x in x_vals:
	tri.append(triaxial_NFW_potential(x,y,z,host))
	# sph.append(spherical_NFW_potential(x,y,z,host_sph))
	tri_gala.append(triaxial.value([x,y,z]))


plt.figure()
plt.plot(x_vals,tri,'c',linewidth=2,label='triaxial')
# plt.plot(x_vals,sph,'b',linewidth=2,label='spherical')
plt.plot(x_vals,tri_gala,'g',linewidth=2,label='gala')
plt.legend()
plt.xscale('log')
plt.show()