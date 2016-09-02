from __future__ import division

import numpy as np

def spherical_coordinates(position):
	r = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
	theta = np.arccos(position[2]/r)
	phi = np.arctan2(position[1],position[0])
	return r, theta, phi

def cartesian_coordinates(r,theta,phi):
	x = r * np.cos(phi) * np.sin(theta)
	y = r * np.sin(phi) * np.sin(theta)
	z = r * np.cos(theta)
	return x, y, z

def rotate(position,axis):
        # THIS ONLY WORKS WITH THE SINGLE AXIS DIRECTION
        _, theta2, phi2 = spherical_coordinates(axis)
        _, theta1, phi1 = spherical_coordinates(np.array([1,0,0]))

        phi = np.pi/2 - phi1
        theta = theta1 - theta2
        psi = phi2 - np.pi/2

        c1 = np.cos(phi)
        c2 = np.cos(theta)
        c3 = np.cos(psi)

        s1 = np.sin(phi)
        s2 = np.sin(theta)
        s3 = np.sin(psi)
        
        R1 = np.array([[c1, -s1, 0],[s1, c1, 0], [0, 0, 1]])
        R2 = np.array([[1, 0, 0],[0, c2, -s2], [0, s2, c2]])
        R3 = np.array([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])

        rot1 = np.dot(R1,position)
        rot2 = np.dot(R2,rot1)
        rot3 = np.dot(R3,rot2)

        return rot3