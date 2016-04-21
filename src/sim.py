from __future__ import division
import numpy as np

from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.timestep.particle_evolution import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.halos.subhalo import *
from sidm_orbit_calculation.src.plotting.make_plots import *

host = HostHalo(M=1.e14*M_sol,potential='point_mass') # changed mass to 1.e13 to match jiang paper! (!!)
gravity = GetGravitationalForce(host)
dt = 1.e11

print 'host halo:'
print 'R200 = ', host.R_200
print 'v200 = ', host.v_200
print 'M = ', host.M

initial_position = [host.R_200,0.]
initial_momentum = [-0.*host.v_200,0.8*host.v_200]

# subhalo = Subhalo(mass_ratio=0.001,initial_position=initial_position,initial_momentum=initial_momentum)
# subhalo class might be unnecessary - combine with position/momentum classes? (!!)

def initiate_particle(position=initial_position,momentum=initial_momentum):
	position = ParticlePosition(initial_position=position,dt=dt)
	momentum = ParticleMomentum(initial_momentum=momentum,initial_force=gravity.calculate_gravitational_force(position=position.current_position),initial_partial_force=gravity.partial_force(position=position.current_position),dt=dt)
	return position,momentum

def step(force=None,dt=None):
	position.update_step(momentum=momentum.current_momentum)
	momentum.update_step(position=position.current_position,force=force)
	return None

def sim(n_orbits=20,dt=1.e14):
	time = 0.
	times = [time]

	while position.current_position[1] < n_orbits*(2*np.pi):
		force = gravity.calculate_gravitational_force(position=position.current_position)
		step(force,dt)
		time+=dt
		times.append(time)
		print 'phi = ', position.current_position[1]/(2*np.pi)
		print 'radius = ', position.current_position[0]

		printing = 0.
		if printing:
			print 'time = ', time
			print 'current position = ', position.current_position
			print 'radius = ', position.current_position[0]
			print 'phi/pi = ', position.current_position[1]/np.pi
			print 'force = ', force
			print 'radial velocity = ', momentum.current_momentum[0]
	
	times = np.array(times)
	positions = np.array(position.position_array)
	momenta = np.array(momentum.momentum_array)
	forces = np.array(gravity.force_array)

	return times, positions, momenta, forces

# FOR NOW MANUALLY DEFINE INITIAL PARAMETERS
position, momentum = initiate_particle(position=initial_position,momentum=initial_momentum)

n_orbits = 20
times, positions, momenta, forces = sim(n_orbits=n_orbits,dt=dt)

printing = 0.
if printing:
	print 'r min = ', positions[:,0].min()/host.R_200
	print 'r max = ', positions[:,0].max()/host.R_200
	print 'r initial = ', positions[:,0][0]/host.R_200
	print 'r final = ', positions[:,0][-1]/host.R_200
	print 'final time = ', time


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
	# plot.angular_position()
	# plot.gravitational_force()
	# plot.radial_velocity()



