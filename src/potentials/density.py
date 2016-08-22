from __future__ import division
import numpy as np

def constant_density(x,y,z,host):
	return 1e-40 # average, just a test

def spherical_NFW_density(x,y,z,host):
	r = np.sqrt(x*x + y*y + z*z)
	R = r/host.R_s
	return host.rho_s/(R*(1+R)**2)

def triaxial_NFW_density(x,y,z,host):
	x2 = x*x
	y2 = y*y
	z2 = z*z
	q2 = host.q*host.q
	s2 = host.s*host.s
	r = np.sqrt(x2 + y2/q2 + z2/s2)
	R = r/host.R_s	
	alpha = 1
	eta = 3
	return host.rho_s / (R**alpha * (1+R)**(eta-alpha))

density_dict = {'point_mass': constant_density, 'spherical_NFW': spherical_NFW_density, 'triaxial_NFW': triaxial_NFW_density}