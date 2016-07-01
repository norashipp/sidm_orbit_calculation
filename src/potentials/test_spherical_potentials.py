import numpy as 
from scipy.special import ellipkinc
from scipy.integrate import quad
from scipy.optimize import brentq

from sidm_orbit_calculation.src.utils.constants import *

'''
Collection of spherically symmetric test potentials
phi = phi(r)
'''

def spherical_harmonic_oscilator(r, A=1, B=1):
    return (A + B*r**2)

def point_mass_potential(r, host):
	ep = 0.
	return -G*host.M/(r+ep)

def isochrone_potential(r, host, b=1):
   return -G*host.M/(b+np.sqrt(b**2+r**2))

def spherical_NFW_potential(r, host):
	return -4*np.pi*G*host.rho_s*host.R_s**2*np.log(1+r/host.R_s)/(r/host.R_s)

def triaxial_NFW_potential(x, y, z, host, s = 1, q = 1):
	R = np.sqrt(x**2+y**2/q**2+z**2/s**2)
	xi = R/host.R_s
	A = 4*np.pi*G*s*q*host.rho_s*host.R_s**2/np.sqrt(1-s**2)

	alpha = 1 # what are alpha and eta again? in this case eta = 3

	phi_out = -A/(2-alpha)*ellipkinc(np.sqrt(1-s**2),np.sqrt((1-q**2)/(1-s**2)))*(1-(xi/(1+xi)**(2-alpha)))
	
	func = lambda lmbda, m: x**2/(m**2*R(x,y,z)**2+lmbda) + y**2/(m**2*R(x,y,z)**2*q**2+lmbda) + z**2/(m**2*R(x,y,z)**2*s**2+lmbda) - 1
	a = -1 # no idea what is reasonable here...
	b = 1
	lmbda = lambda m: brentq(func,a,b,args=(m)) # is this how args works
	# lmbda = lambda m: # not sure how to find the root of this equation...
	integrand = lambda m: m**(1-alpha)/((1+m*xi)**(3-alpha))*ellipkinc(np.sqrt((1-s**2)/(1+lmbda(m)/(m**2)))
	integral = quad(integrand,0,1)
	phi_in = -A * xi ** (2 - alpha) * integral 
	
	return phi_out + phi_in

potential_dict = {'spherical_harmonic':spherical_harmonic_oscilator,'point_mass':point_mass_potential,'isochrone':isochrone_potential,'spherical_NFW':spherical_NFW_potential}
