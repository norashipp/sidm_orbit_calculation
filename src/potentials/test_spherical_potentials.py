import numpy as np

from sidm_orbit_calculation.src.utils.constants import *

'''
Collection of spherically symmetric test potentials
phi = phi(r)
'''

def spherical_harmonic_oscilator(r, A=1, B=1) :
    return (A + B*r**2)

def point_mass_potential(r, host) :
	ep = 0.
	return -G*host.M/(r+ep)

def isochrone_potential(r, host, b=1) :
   return -G*host.M/(b+np.sqrt(b**2+r**2))

def spherical_NFW_potential(r, host):
	return -4*np.pi*G*host.rho_s*host.R_s**2*np.log(1+r/host.R_s)/(r/host.R_s)

potential_dict = {'spherical_harmonic':spherical_harmonic_oscilator,'point_mass':point_mass_potential,'isochrone':isochrone_potential,'spherical_NFW':spherical_NFW_potential}
