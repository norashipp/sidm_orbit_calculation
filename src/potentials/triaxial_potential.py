import numpy as np
from scipy.special import ellipkinc
from scipy.integrate import quad
from scipy.optimize import brentq

from sidm_orbit_calculation.src.utils.constants import *


# IN PROGRESS

def calculate_lambda(m,x,y,z,host):
	R = radius(x,y,z,host)
	func = lambda L: x**2/(m**2*R**2+L) + y**2/(m**2*R**2*host.q**2+L) + z**2/(m**2*R**2*host.s**2+L) - 1
	vals = np.linspace(-R**2,R**2)
	root = -1*np.inf
	for i in range(len(vals)-1):
		ratio = func(vals[i])/func(vals[i+1])
		if ratio < 0 and i > root:
			root = ratio
			a = vals[i]
			b = vals[i+1]
	lmbda = brentq(func,a,b)
	return lmbda

def integrand(m,x,y,z,xi,alpha,host):
	lmbda = calculate_lambda(m,x,y,z,host)
	return m**(1-alpha)/((1+m*xi)**(3-alpha))*ellipkinc(np.sqrt((1-host.s**2)/(1+lmbda/(m**2*radius(x,y,z,host)**2))),np.sqrt((1-host.q**2)/(1-host.s**2)))
	
def inner_potential(x,y,z,xi,A,alpha,host):
	integral = quad(integrand,0,1,args=(x,y,z,xi,alpha,host))[0]
	return -A * xi ** (2 - alpha) * integral 

def outer_potential(xi,A,alpha,host):
	return -A/(2-alpha)*ellipkinc(np.sqrt(1-host.s**2),np.sqrt((1-host.q**2)/(1-host.s**2)))*(1-(xi/(1+xi)**(2-alpha)))
	
def radius(x,y,z,host):
	return np.sqrt(x**2+y**2/host.q**2+z**2/host.s**2)

def triaxial_NFW_potential(position,host):
	x, y, z = position
	xi = radius(x,y,z,host)/host.R_s
	A = 4*np.pi*G*host.s*host.q*host.rho_s*host.R_s**2/np.sqrt(1-host.s**2)

	alpha = 1 # what are alpha and eta again? in this case eta = 3

	inner = inner_potential(x,y,z,xi,A,alpha,host)
	outer = outer_potential(xi,A,alpha,host)
	# phi = inner_potential(x,y,z,xi,A,alpha,host) + outer_potential(xi,A,alpha,host)
	phi = inner + outer
	# print inner
	# print outer
	# print phi
	return phi