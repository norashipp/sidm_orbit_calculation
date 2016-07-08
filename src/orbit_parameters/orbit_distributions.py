'''Set up a distribution of orbit parameters'''

# for now only allowing M_h = 1e12, 1e13, 1e14. later interpolate in that dimension as well.

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
	return gaussian(x = x, sigma = sigma, mu = mu)*lorentz(x = z, gamma = gamma)

def total_velocity(sigma, gamma, mu):
	# x = V/V200
	res = scipy.integrate.quad(integrand, -np.inf, np.inf, agrs = (x, sigma, gamma, mu))[0]
	return res

def radial_ratio(B):
	return A * (np.exp(B * V_r / V) - 1)
	# A is a normalization constant

##########################################

def jiang14():
	B, sigma, gamma, mu = fitting_paramters(M_h = M_h, M_s = M_s)

