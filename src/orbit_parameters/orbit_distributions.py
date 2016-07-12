import numpy as np
from scipy.integrate import quad
from numpy import random
from scipy.interpolate import UnivariateSpline

'''Set up a distribution of orbit parameters'''

# For now only allowing M_h = 1e12, 1e13, 1e14. later interpolate in that dimension as well.

def read_data(fname):
	f = open(fname, 'r')
	lines = f.readlines()[1:]
	f.close()
	x = []
	y = []
	for line in lines:
	    line = line.strip()
	    line = line.split(',')
	    x.append(float(line[0]))
	    y.append(float(line[1]))
	return UnivariateSpline(x, y, s=0, k=1)

def fitting_parameters(M_h, M_s):
	# use data thief to get fitting parameters as a funciton of M_h and M_s/M_h
	x = np.log10(M_s / M_h)
	
	params = ['B', 'sigma', 'gamma', 'mu']
	param_vals = []
	
	masses = ['1e12', '1e13', '1e14']
	mass_vals = [1e12, 1e13, 1e14]
	m = masses[np.argmin(np.abs(M_h - mass_vals))]

	file_name = 'data_thief_' + m + '_'

	for p in params:
		fname = file_name + p + '.txt'
		sp = read_data(fname = fname)
		param_vals.append(sp(x))

	return param_vals

##########################################

def lorentz(x, gamma):
	return gamma / (np.pi * (x ** 2 + gamma ** 2))

def gaussian(x, sigma, mu):
	return 1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

def integrand(y, x, sigma, gamma, mu):
	# integrating over y
	z = x - y
	return gaussian(x, sigma, mu) * lorentz(x,gamma)

def total_velocity(x, sigma, gamma, mu):
	# x = V/V200
	res = scipy.integrate.quad(integrand, -np.inf, np.inf, agrs = (x, sigma, gamma, mu))[0]
	return res

def radial_velocity(x, A, B):
	# A is a normalization constant
	# x = Vr/V
	return A * (np.exp(B * x) - 1)


##########################################

def normalize(B):
	func = lambda x: np.exp(B*x)-1
	res = scipy.integrate.quad(func,0,1)[0]
	A = 1/res
	return A

def invert_radial(A,B):
	x = np.linspace(0, 1, 100)
	y = radial_velocity(x, A,  B)
	sp = UnivariateSpline(y, x, s=0, k=1)
	return sp

def invert_total(sigma, gamma, mu):
	x = np.linspace(0, 3.0, 100) # bounds taken from jiang fig 6
	y = total_velocity(x, simga, gamma, mu)
	sp = UnivariateSpline(y, x, s=0, k=1)
	return sp

##########################################

def get_inverse_cdf(M_h, M_s):
	# x1, x2 pair of random numbers [0,1]
	B, sigma, gamma, mu = fitting_paramters(M_h, M_s)
	A = normalize(B)
	total_inverse_cdf = invert_total(sigma, gamma, mu)
	radial_inverse_cdf = invert_radial(A, B)
	return total_inverse_cdf, radial_inverse_cdf

def calculate_velocities(total_ratio, radial_ratio, host):
	v_r = radial_ratio * total_radtio * host.v_200
	v_theta = total_ratio * np.sqrt(1 - radial_ratio ** 2)
	return v_r, v_theta

def get_velocities(subhalo):
	# is there a way to reduce computation here? need to recalculate fitting parameters for each subhalo mass, though
	x_ total, x_radial = np.random.uniform(0,1,2)
	total_inverse_cdf, radial_inverse_cdf = get_inverse_cdf(subhalo.host.M, subhalo.M)
	total_ratio = total_inverse_cdf(x_total)
	radial_ratio = radial_inverse_cdf(x_radial)
	v_r, v_theta = calculate_velocities(vt_ratio, vr_ratio)
	return v_r, v_theta














