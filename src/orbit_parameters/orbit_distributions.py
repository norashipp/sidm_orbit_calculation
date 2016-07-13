import numpy as np
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline
from numpy.random import choice
from numpy.random import uniform

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
	mass_vals = np.array([1e12, 1e13, 1e14])
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
	integrand(0, x, sigma, gamma, mu)
	res = quad(integrand, -np.inf, np.inf, args = (x, sigma, gamma, mu))[0]
	return res

def radial_velocity(x, A, B):
	# A is a normalization constant
	# x = Vr/V
	return A * (np.exp(B * x) - 1)


##########################################

def normalize(B):
	func = lambda x: np.exp(B*x)-1
	res = quad(func,0,1)[0]
	A = 1/res
	return A

def invert_radial(A,B):
	xx = np.linspace(0, 1, 100)
	yy = []
	for x in xx:
		y = radial_velocity(x, A,  B)
		yy.append(y)
	sp = UnivariateSpline(yy, xx, s=0, k=1)
	return sp

def invert_total(sigma, gamma, mu):
	xx = np.linspace(0, 3.0, 100) # bounds taken from jiang fig 6
	yy = []
	for x in xx:
		y = total_velocity(x, sigma, gamma, mu)
		yy.append(y)
	sp = UnivariateSpline(yy, xx, s=0, k=1)
	return sp

##########################################

def get_inverse_cdf(M_h, M_s):
	# x1, x2 pair of random numbers [0,1]
	B, sigma, gamma, mu = fitting_parameters(M_h, M_s)
	A = normalize(B)
	total_inverse_cdf = invert_total(sigma, gamma, mu)
	radial_inverse_cdf = invert_radial(A, B)
	return total_inverse_cdf, radial_inverse_cdf

# def calculate_velocities(total_ratio, radial_ratio):
# 	v_r = radial_ratio * total_radtio
# 	v_theta = total_ratio * np.sqrt(1 - radial_ratio ** 2)
# 	return v_r, v_theta

def cartesian_velocities(total_ratio, radial_ratio):
	vx = radial_ratio * total_ratio
	vy = total_ratio * np.sqrt(1 - radial_ratio ** 2)
	vz = 0
	return np.array([vx, vy, vz])

def rotate_orbit(velocity):
	u = np.random.uniform(0,1)
	v = np.random.uniform(0,1)

	phi = 2*np.pi*u
	theta = np.arccos(2*v - 1)
	
	# cos_theta = np.random.uniform(-1,1)
	# phi = np.random.uniform(0,2*np.pi)

	# sin_theta = np.sqrt(1-cos_theta**2)
	north_south = choice([-1,1])
	if north_south == -1: theta+=np.pi
	# sin_theta *= north_south

	# rotate around x by theta
	# rotate around z by phi
	# cos theta uniformly distributed
	# phi uniformly distributed
	# randomly select between northern and southern hemisphere

	Rx = np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
	Rz = np.array([[np.cos(phi), -np.sin(phi), 0], [np.sin(phi), np.cos(phi), 0], [0, 0, 1]])

	velocity_rotated = np.dot(Rx,velocity)
	velocity_rotated = np.dot(Rz, velocity_rotated)

	position = np.array([1,0,0])
	position_rotated = np.dot(Rx, position)
	position_rotated = np.dot(Rz, position_rotated)

	return position_rotated, velocity_rotated

def initial_conditions(subhalo):
	# is there a way to reduce computation here? need to recalculate fitting parameters for each subhalo mass, though
	
	total_inverse_cdf, radial_inverse_cdf = get_inverse_cdf(subhalo.host.M, subhalo.M)
	
	x_total, x_radial = uniform(0,1,2)

	total_ratio = total_inverse_cdf(x_total)
	radial_ratio = radial_inverse_cdf(x_radial)
	
	velocity = cartesian_velocities(total_ratio, radial_ratio)

	initial_position, initial_velocity = rotate_orbit(velocity)

	return initial_position*subhalo.host.R_200, initial_velocity*subhalo.host.v_200














