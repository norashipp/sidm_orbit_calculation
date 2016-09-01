from __future__ import division

import numpy as np

def spherical_coordinates(position):
	r = np.sqrt(position[0]**2+position[1]**2+position[2]**2)
	theta = np.arccos(position[2]/r)
	phi = np.arctan2(position[1],position[0])
	return r, theta, phi

def cartesian_coordinates(r,theta,phi):
	x = r * np.cos(phi) * np.sin(theta)
	y = r * np.sin(phi) * np.sin(theta)
	z = r * np.cos(theta)
	return x, y, z