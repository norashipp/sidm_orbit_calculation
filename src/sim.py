from __future__ import division
import numpy as np
import cPickle

from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.timestep.particle_evolution import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.halos.subhalo import *
from sidm_orbit_calculation.src.plotting.make_plots import *

host = HostHalo(M=1.e13*M_sol,potential='point_mass') # changed mass to 1.e13 to match jiang paper! (!!)
gravity = GetGravitationalForce(host)
dt = 1.e12

print 'host halo:'
print 'R200 = ', host.R_200
print 'v200 = ', host.v_200
print 'M = ', host.M

r0 = host.R_200
phi0 = 0.
initial_position = [r0*np.cos(phi0),r0*np.sin(phi0)]
initial_position = [r0,phi0]
initial_momentum = [0.,1.e6]

def initiate_particle(position=initial_position,momentum=initial_momentum):
	position = ParticlePosition(initial_position=position,dt=dt)
	momentum = ParticleMomentum(initial_momentum=momentum,dt=dt)
	# momentum = ParticleMomentum(initial_momentum=momentum,initial_force=gravity.calculate_gravitational_force(position=position.current_position),initial_partial_force=gravity.partial_force(position=position.current_position),dt=dt)
	return position,momentum

def step(dt=None):
	position.update_step(momentum=momentum.current_momentum)
	momentum.update_step(position=position.current_position,force=gravity.calculate_gravitational_force(position=position.current_position))
	return None

def sim(n_orbits=20,dt=1.e14):
	time = 0.
	times = [time]

	orbit_time = time

	while position.current_position[1] < n_orbits*(2*np.pi) and time < dt*1.e6:
		step(dt)
		time+=dt
		times.append(time)
		
		# print 'phi = ', position.current_position[1]
		# print 'phi fraction = ', position.current_position[1]/(2*np.pi)
		# print 'radius = ', position.current_position[0]
		# print 'x = ', position.current_position[0]*np.cos(position.current_position[1])
		# print 'y = ', position.current_position[0]*np.sin(position.current_position[1])
		# assert 0

		printing = 1.
		if printing:
			print 'time = ', time
			print 'position =  %10.5g %10.5g' % (position.current_position[0],position.current_position[1])
			print 'force = %10.5g' % gravity.calculate_gravitational_force(position=position.current_position)[0]
			# assert 0
	
	times = np.array(times)
	positions = np.array(position.position_array)
	momenta = np.array(momentum.momentum_array)
	forces = np.array(gravity.force_array)

	return times, positions, momenta, forces

def write_output(data):
	f = open('data/pickle.dat','wb')
	cPickle.dump(data,f)
	f.close()
	return None

# FOR NOW MANUALLY DEFINE INITIAL PARAMETERS
position, momentum = initiate_particle(position=initial_position,momentum=initial_momentum)

n_orbits = 2
times, positions, momenta, forces = sim(n_orbits=n_orbits,dt=dt)
write_output([times,positions,momenta,forces,host])

printing = 1.
if printing:
	print 'r min = ', positions[:,0].min()/host.R_200
	print 'r max = ', positions[:,0].max()/host.R_200
	print 'r initial = ', positions[:,0][0]/host.R_200
	print 'r final = ', positions[:,0][-1]/host.R_200
	print 'final time = ', times.max()


dt_test = 0
if dt_test:
	f = open('dt_test.txt','a')
	f.write(str(dt)+','+str(positions[:,0][-1]/host.R_200)+'\n')
	f.close()

dx_test = 0
if dx_test:
	f = open('dx_test.txt','a')
	f.write(str(dxdiv)+','+str(positions[:,0][-1]/host.R_200)+'\n')
	f.close()

plotting = 1
if plotting:
	plot = Plotting(times=times,positions=positions,momenta=momenta,forces=forces,host=host)
	plot.orbit()
	# plot.orbit_color()
	plot.radial_position()
	plot.angular_position()
	plot.gravitational_force()
	# plot.radial_velocity()



