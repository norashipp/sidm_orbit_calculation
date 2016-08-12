import numpy as np
from scipy.special import ellipkinc
from scipy.integrate import quad
from scipy.optimize import brentq
import memory_profiler

from sidm_orbit_calculation.src.utils.constants import *

def calculate_lambda(m, x, y, z, host):
	R = radius(x, y, z, host)
	x2, y2, z2 = x*x, y*y, z*z
	mr2, mrq2, mrs2 = m*m*R*R, m*m*R*R*host.q*host.q, m*m*R*R*host.s*host.s
	func = lambda L: x2/(mr2+L) + y2/(mrq2+L) + z2/(mrs2+L) - 1
	root = -1*np.inf

	a, b = calculate_range(func, R, root)

	lmbda = brentq(func, a, b)
	return lmbda

VALS_LEN = 50
INDICES = np.arange(VALS_LEN - 1, dtype=int)
EPS = 1e6
count = 0
prev_R, prev_vals = None, None

def almost_eq(x, y):
	return np.abs(x - y) / x < x / EPS

def calculate_range(func, R, root):
	global prev_R, count, prev_vals
	if prev_R == R:
		vals = prev_vals
	else:
		prevR = R
		vals = np.linspace(-R**2, R**2, VALS_LEN)
		prev_vals = vals

	fs = func(vals)

	ratios = fs[:-1]/fs[1:]
	i = INDICES[ratios < 0][-1]
	a, b = vals[i], vals[i+1]

	return a, b

def integrand(m,x,y,z,xi,alpha,host):
	lmbda = calculate_lambda(m, x, y, z, host)
	s2, m2, q2 = host.s*host.s, m*m, host.q*host.q
	rad = radius(x, y, z, host)
	return m**(1 - alpha) / ((1 + m*xi)**(3 - alpha)) * ellipkinc(np.sqrt((1 - s2)/(1 + lmbda/(m2 * rad*rad))), np.sqrt((1 - q2)/(1 - s2)))
	
def inner_potential(x,y,z,xi,A,alpha,host):
	integral = quad(integrand, 0, 1, args=(x, y, z, xi, alpha, host))[0]
	return -A * xi**(2 - alpha) * integral 

def outer_potential(xi,A,alpha,host):
	return -A/(2 - alpha) * ellipkinc(np.sqrt(1 - host.s*host.s), np.sqrt((1 - host.q*host.q)/(1 - host.s*host.s))) * (1 - ((xi/(1 + xi))**(2 - alpha)))
	
def radius(x,y,z,host):
	return np.sqrt(x*x + y*y / (host.q*host.q) + z*z / (host.s*host.s))

def triaxial_NFW_potential(x, y, z, host): # not working within 1e18 m - that's 0.0001 times the virial radius, is this ok? smoothing?
	# x, y, z = position
	xi = radius(x, y, z, host) / host.R_s
	A = 4 * np.pi * G * host.s * host.q * host.rho_s * host.R_s**2 / np.sqrt(1 - host.s**2)

	alpha = 1 # eta = 3

	inner = inner_potential(x,y,z,xi,A,alpha,host)
	outer = outer_potential(xi,A,alpha,host)
	# phi = inner_potential(x,y,z,xi,A,alpha,host) + outer_potential(xi,A,alpha,host)
	phi = inner + outer
	# print inner
	# print outer
	# print phi
	return phi
