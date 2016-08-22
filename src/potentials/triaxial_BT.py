from __future__ import division

import numpy as np
from scipy.integrate import quad

from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.utils.geometry import *

alpha = 1

def F1(u):
	# integrand = lambda t: (t/(t+1))**(2-alpha)
	# integral = quad(integrand,0,u)[0]

	# return 1/(alpha-2) * (1-1/u * integral)

	return -np.log(1+u)/u

def F2(u):
	# integrand1 = lambda t: (t/(t+1))**(2-alpha)
	# integral1 = quad(integrand1,0,u)[0]

	# integrand2 = lambda t: t**(4-alpha) / (t+1)**(2-alpha) 
	# integral2 = quad(integrand2,0,u)[0]

	# return 1/(alpha-2) * (-2/3 + 1/u*integral1 - 1/u**3*integral2)

	return -1./3 + (2*u*u - 3*u + 6)/(6*u*u) + (1/u - u**-3)*np.log(1+u)

def F3(u):
	# integrand = lambda t: t**(4-alpha) / (t+1)**(3-alpha)
	# integral = quad(integrand,0,u)[0]

	# return 1/u**3*integral

	return (u*u - 3*u - 6)/(2*u*u*(1+u)) + 3*u**-3*np.log(1+u)

def triaxial_NFW_potential(x, y, z, host):
	# r, theta, phi = spherical_coordinates([x,y,z])
	r = np.sqrt(x**2 + y**2 + z**2)

	# eb = np.sqrt(1 - host.q**2)
	# ec = np.sqrt(1 - host.s**2)

	# a = 1
	# b = host.q
	# c = host.s

	# eb2 = 1-(b/a)**2
	# ec2 = 1-(c/a)**2

	eb2 = 1 - host.q**2
	ec2 = 1 - host.s**2

	# M = host.mass_function(host=host, a=0, b=host.R_s)
	# v_c = np.sqrt(G*M/host.R_s)

	# M = 2.10670171738e+13
	# v_c = 0.810346103989
	# R = 0.144898
	# v_c = np.sqrt(G*M/host.R_s)
	# a = 1
	# b = 0.58406
	# c = 0.44725

	# phi0 = v_c*v_c / (np.log(2) - 0.5 + (np.log(2)-0.75)*eb2 + (np.log(2)-0.75)*ec2)

	u = r/host.R_s
	# u = r/R

	C = 4*np.pi*G * host.rho_s * host.R_s**2
	
	costh2 = z*z / (r*r)
	sinth2 = 1 - costh2 
	sinph2 = y*y / (x*x + y*y)

	pot = C*(F1(u) + (eb2 + ec2)/2 * F2(u) + (eb2 * sinth2 * sinph2 + ec2 * costh2)/2 * F3(u))
	
	# pot = phi0 * (F1(u) + (eb2+ec2)/2*F2(u) + (eb2*sinth2*sinph2 + ec2*costh2)/2 * F3(u))

	return pot
