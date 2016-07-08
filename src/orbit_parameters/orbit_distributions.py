'''Set up a distribution of orbit parameters'''

def fitting_parameters(M_h,M_s):
	x = M_s/M_h
	# use data thief to get fitting parameters as a funciton of M_h and M_s/M_h
	return B, sigma, gamma, mu

def lorentz():
	return gamma/(np.pi*(x**2+gamma**2))

def gaussian():
	return 1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


def total_velocity():
	integrand = gaussian()*lorentz()
	res = scipy.integrate.quad(integrand,-np.inf,np.inf)[0]
	return res

def radial_ratio():
	return A * (np.exp(B * V_r / V) - 1)

def jiang14():
	