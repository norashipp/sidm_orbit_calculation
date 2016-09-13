from __future__ import division
import numpy as np
import time

from sidm_orbit_calculation.src.potentials.test_spherical_potentials import *

def euler(subhalo,dt,sigma):
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
        
def leapfrog(subhalo,dt,sigma):
	x0 = subhalo.position[:]
	p0 = subhalo.momentum[:]

	_, f0 = subhalo.gravity.calculate_gravitational_force(position=np.copy(x0))

	x1 = x0 + p0 * dt + 0.5 * f0 * dt ** 2

	_, f1 = subhalo.gravity.calculate_gravitational_force(position=np.copy(x1))

	p1 = p0 + 0.5 * (f0 + f1) * dt
	
	'''
	print 'leapfrog'
	print 'x0: %.3g, %.3g, %.3g' %(x0[0],x0[1],x0[2])
	print 'p0: %.3g, %.3g, %.3g' %(p0[0],p0[1],p0[2])
	print 'r = ', np.sqrt(x0[0]**2+x0[1]**2+x0[2]**2)
	print 'f0: %.3g, %.3g, %.3g' %(f0[0],f0[1],f0[2])
	print 'x1: %.3g, %.3g, %.3g' %(x1[0],x1[1],x1[2])
	print 'p1: %.3g, %.3g, %.3g' %(p1[0],p1[1],p1[2])
	# print '%.3g, %.3g, %.3g' %(x1[0]-x0[0]-1e14,x1[1]-x0[1]-1e14,x1[2]-x0[2]-1e14)
	print 
	time.sleep(2)
	'''
	
	subhalo.position = x1
	subhalo.momentum = p1

	subhalo.gravity.force = f1
	subhalo.drag.force = np.zeros_like(subhalo.position)

	update_arrays(subhalo=subhalo)

def dissipative(subhalo,dt,sigma):
	lmbda = 0.5

	x0 = subhalo.position[:]
	p0 = subhalo.momentum[:]
	
	_, fg0 = subhalo.gravity.calculate_gravitational_force(position=np.copy(x0))
	
	fd0 = subhalo.drag.calculate_drag_force(position=np.copy(x0),momentum=np.copy(p0),subhalo_mass=subhalo.M)
	# print 'rho0: %.3g' %subhalo.host.rho 
	
	f0 = fg0 + fd0

	x1 = x0 + dt * p0 + 0.5 * dt ** 2 * f0

	p_half = p0 + lmbda * dt * f0

	_, fg1 = subhalo.gravity.calculate_gravitational_force(position=np.copy(x1))
	fd1 = subhalo.drag.calculate_drag_force(position=np.copy(x0),momentum=np.copy(p0),subhalo_mass=subhalo.M)
	f1 = fg1 + fd1

	p1 = p0 + 0.5 * dt * (f0 + f1)

	'''
	print 'dissipative'
	# print 'x0: %.3g, %.3g, %.3g' %(x0[0],x0[1],x0[2])
	print 'p0: %.3g, %.3g, %.3g' %(p0[0],p0[1],p0[2])
	print 'f0: %.3g, %.3g, %.3g' %(f0[0],f0[1],f0[2])
	print 'fd0: %.3g, %.3g, %.3g' %(fd0[0], fd0[1], fd0[2])
	print 'fg0: %.3g, %.3g, %.3g' %(fg0[0], fg0[1], fg0[2])
	# print 'x1: %.3g, %.3g, %.3g' %(x1[0],x1[1],x1[2])
	# print 'p1: %.3g, %.3g, %.3g' %(p1[0],p1[1],p1[2])
	print x1-x0
	print p1-p0
	print 
	time.sleep(2)
	'''

	subhalo.position = x1
	subhalo.momentum = p1
	subhalo.drag.force = fd1
	subhalo.gravity.force = fg1

	update_arrays(subhalo=subhalo)

integrator_dict = {'euler':euler,'leapfrog':leapfrog,'dissipative':dissipative}

#################################

def update_arrays(subhalo):
	# if not subhalo.count % 500:
	subhalo.position_array.append((subhalo.position[0],subhalo.position[1],subhalo.position[2]))
	subhalo.momentum_array.append((subhalo.momentum[0],subhalo.momentum[1],subhalo.momentum[2]))
	subhalo.gravity.force_array.append(subhalo.gravity.force)
	subhalo.drag.force_array.append(subhalo.drag.force)
	subhalo.host.density_array.append(subhalo.host.rho)
