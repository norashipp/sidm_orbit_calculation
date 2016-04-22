import numpy as np

from sidm_orbit_calculation.src.utils.constants import *

'''
Collection of spherically symmetric test potentials
phi = phi(r)
'''

def spherical_harmonic_oscilator(r=None, A=1, B=1) :
    return -(A + B*r**2)

def point_mass_potential(r=None, host=None) :
	ep = 1e-2*host.R_200
	return G*host.M/(r+ep)

def isochrone_potential(r=None, b=1, host=None) :
   return G*host.M/(b+np.sqrt(b**2+r**2))

def spherical_NFW_potential(r=None, host=None):
	return 4*np.pi*G*host.scale_density*host.scale_radius**2*np.log(1+r/host.scale_radius)/(r/host.scale_radius)

potential_dict = {'spherical_harmonic':spherical_harmonic_oscilator,'point_mass':point_mass_potential,'isochrone':isochrone_potential,'spherical_NFW':spherical_NFW_potential}