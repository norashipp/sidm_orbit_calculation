from __future__ import division
import numpy as np
import time

def euler(subhalo,dt):
	x0 = subhalo.position[:]
	p0 = subhalo.momentum[:]
	f0 = update_gravity(gravity=subhalo.gravity,position=x0)

	x1 = x0 + p0 * dt
	p1 = p0 + f0 * dt

	f1 = update_gravity(gravity=subhalo.gravity,position=x1) # so udpated force array corresponds to udpated position array

	subhalo.position = x1
	subhalo.momentum = p1

	subhalo.drag.force = 0.
	
	update_arrays(subhalo=subhalo)
        
def leapfrog(subhalo,dt):
	x0 = subhalo.position[:]
	p0 = subhalo.momentum[:]
	
	f0 = update_gravity(gravity=subhalo.gravity,position=x0)

	x1 = x0 + p0 * dt + 0.5 * f0 * dt ** 2

	f1 = update_gravity(gravity=subhalo.gravity,position=x1)

	p1 = p0 + 0.5 * (f0 + f1) * dt
	
	subhalo.position = x1
	subhalo.momentum = p1

	subhalo.drag.force = np.zeros_like(subhalo.position)
	subhalo.host.update_density(position=subhalo.position)

	update_arrays(subhalo=subhalo)

def dissipative(subhalo,dt):
	lmbda = 0.5

	x0 = subhalo.position[:]
	p0 = subhalo.momentum[:]
	fg0 = update_gravity(gravity=subhalo.gravity,position=x0)
	fd0 = subhalo.drag.calculate_drag_force(position=subhalo.position,momentum=subhalo.momentum)
	f0 = fg0 + fd0

	x1 = x0 + dt * p0 + 0.5 * dt ** 2 * f0

	p_half = p0 + lmbda * dt * f0

	fg1 = update_gravity(gravity=subhalo.gravity,position=x1)
	fd1 = subhalo.drag.calculate_drag_force(position=subhalo.position,momentum=subhalo.momentum) # should this be in update force function?
	# print 'drag force = ', fd1
	f1 = fg1 + fd1

	p1 = p0 + 0.5 * dt * (f0 + f1)

	subhalo.position = x1
	subhalo.momentum = p1
	# subhalo.gravity.force = fg1 # figure out this updating - updated in get grav
	subhalo.drag.force = fd1

	update_arrays(subhalo=subhalo)

integrator_dict = {'euler':euler,'leapfrog':leapfrog,'dissipative':dissipative}

#################################

def update_arrays(subhalo):
	if not subhalo.count % 500:
		subhalo.position_array.append((subhalo.position[0],subhalo.position[1],subhalo.position[2]))
		subhalo.momentum_array.append((subhalo.momentum[0],subhalo.momentum[1],subhalo.momentum[2]))
		subhalo.gravity.force_array.append(subhalo.gravity.force)
		subhalo.drag.force_array.append(subhalo.drag.force)
		subhalo.host.density_array.append(subhalo.host.rho)

# UPDATE TO 3D
def update_gravity(gravity,position):
    gravity.calculate_gravitational_force(position=position)
    r, theta, phi = spherical_coordinates(position=position)
    force = gravity.vector(phi=phi,theta=theta)
    # print 'gravitational force = ', force
    return force

def spherical_coordinates(position):
	r = np.sqrt(position[0]**2+position[1]**2+position[2]**2)
	theta = np.arccos(position[2]/r)
	phi = np.arctan2(position[1],position[0])
	# print 'x, y , z = ', position
	# print 'r, theta, phi = ', r, theta, phi
	return r, theta, phi

