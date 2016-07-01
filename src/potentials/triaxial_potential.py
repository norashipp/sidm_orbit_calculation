import numpy as 
from scipy.special import ellipkinc
from scipy.integrate import quad
from scipy.optimize import brentq

from sidm_orbit_calculation.src.utils.constants import *


# IN PROGRESS

def innner_potential():
	func = lambda lmbda, m: x**2/(m**2*R(x,y,z)**2+lmbda) + y**2/(m**2*R(x,y,z)**2*q**2+lmbda) + z**2/(m**2*R(x,y,z)**2*s**2+lmbda) - 1
	a = -1 # no idea what is reasonable here...
	b = 1
	lmbda = lambda m: brentq(func,a,b,args=(m)) # is this how args works
	# lmbda = lambda m: # not sure how to find the root of this equation...
	integrand = lambda m: m**(1-alpha)/((1+m*xi)**(3-alpha))*ellipkinc(np.sqrt((1-s**2)/(1+lmbda(m)/(m**2)))
	integral = quad(integrand,0,1)
	return -A * xi ** (2 - alpha) * integral 

def outer_potential():
	return -A/(2-alpha)*ellipkinc(np.sqrt(1-s**2),np.sqrt((1-q**2)/(1-s**2)))*(1-(xi/(1+xi)**(2-alpha)))
	
def triaxial_NFW_potential(x, y, z, host, s = 1, q = 1):
	R = np.sqrt(x**2+y**2/q**2+z**2/s**2)
	xi = R/host.R_s
	A = 4*np.pi*G*s*q*host.rho_s*host.R_s**2/np.sqrt(1-s**2)

	alpha = 1 # what are alpha and eta again? in this case eta = 3

	phi = inner_potential(R,xi,A,alpha) + outer_potential(xi,A,alpha)
	return phi