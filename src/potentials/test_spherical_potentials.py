import time
import numpy as np

from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.potentials.triaxial_potential import triaxial_NFW_potential

'''
Collection of spherically symmetric test potentials
phi = phi(r)
'''

def spherical_harmonic_oscilator(x, y, z, A=1, B=1):
    r = np.sqrt(x**2 + y**2 + z**2)
    # r = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
    return (A + B*r**2)

def point_mass_potential(x, y, z, host):
	r = np.sqrt(x**2 + y**2 + z**2)
	# r = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
	ep = 0.
	return -G*host.M/(r+ep)

def isochrone_potential(x, y, z, host, b=1):
	r = np.sqrt(x**2 + y**2 + z**2)
	# r = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
	return -G*host.M/(b+np.sqrt(b**2+r**2))

def spherical_NFW_potential(x, y, z, host):
	r = np.sqrt(x**2 + y**2 + z**2)
	# r = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
	ep = 0.05 # ok with max v0 = -5e6, dt = 1e2
	return -4*np.pi*G*host.rho_s*host.R_s**2*np.log(1+r/host.R_s)/(r/host.R_s+ep)

potential_dict = {'spherical_harmonic':spherical_harmonic_oscilator,'point_mass':point_mass_potential,'isochrone':isochrone_potential,'spherical_NFW':spherical_NFW_potential,'triaxial_NFW':triaxial_NFW_potential}
