from __future__ import division
import numpy as np

def euler(subhalo=None,dt=None):
	x0 = subhalo.position[:]
	p0 = subhalo.momentum[:]
	f0 = update_force(gravity=subhalo.gravity,position=x0)

	x1 = x0 + p0 * dt
	p1 = p0 + f0 * dt

	f1 = update_force(gravity=subhalo.gravity,position=x1) # so udpated force array corresponds to udpated position array

	subhalo.position = x1
	subhalo.momentum = p1
	
	update_arrays(subhalo=subhalo)
        
def leapfrog(subhalo=None,dt=None):
	x0 = subhalo.position[:]
	p0 = subhalo.momentum[:]
	f0 = update_force(gravity=subhalo.gravity,position=x0)

	x1 = x0 + p0 * dt + 0.5 * f0 * dt ** 2

	f1 = update_force(gravity=subhalo.gravity,position=x1)

	p1 = p0 + 0.5 * (f0 + f1) * dt
	
	subhalo.position = x1
	subhalo.momentum = p1

	update_arrays(subhalo=subhalo)

def dissipative(subhalo=None,dt=None):
	lmbda = 0.5

	x0 = subhalo.position[:]
	p0 = subhalo.momentum[:]
	fg0 = update_force(gravity=subhalo.gravity,position=x0)
	fd0 = subhalo.drag.calculate_drag_force(position=subhalo.position,momentum=subhalo.momentum)
	f0 = fg0 + fd0

	x1 = x0 + dt * p0 + 0.5 * dt ** 2 * f0

	p_half = p0 + lmbda * dt * f0

	fg1 = update_force(gravity=subhalo.gravity,position=x1)
	fd1 = subhalo.drag.calculate_drag_force(position=subhalo.position,momentum=subhalo.momentum)
	f1 = fg1 + fd1

	p1 = p0 + 0.5 * dt * (f0 + f1)

	subhalo.position = x1
	subhalo.momentum = p1
	# subhalo.gravity.force = fg1 # figure out this updating
	subhalo.drag.force = fd1

	update_arrays(subhalo=subhalo)

integrator_dict = {'euler':euler,'leapfrog':leapfrog,'dissipative':dissipative}

#################################

def update_arrays(subhalo=None):
	subhalo.position_array.append((subhalo.position[0],subhalo.position[1]))
	subhalo.momentum_array.append((subhalo.momentum[0],subhalo.momentum[1]))
	subhalo.gravity.force_array.append(subhalo.gravity.force)
	subhalo.drag.force_array.append(subhalo.drag.force)

def update_force(gravity=None,position=None):
    gravity.calculate_gravitational_force(position=position)
    phi = calculate_angle(position=position)
    force = gravity.vector(phi=phi)
    return force

def calculate_angle(position=None):
    return np.arctan2(position[1],position[0])

