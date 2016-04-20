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

print 'host halo:'
print 'R200 = ', host.R_200
print 'v200 = ', host.v_200
print 'M = ', host.M

# orbits = OrbitParameters()

# momenta = orbits.momentum()
# positions = orbits.position()
# for i in n_subhalos:
# 	subhalos.append(Subhalos(mass_ratio=,initial_position=positions[i],initial_momentum=momenta[i]))

# FOR NOW MANUALLY DEFINE INITIAL PARAMETERS

initial_position = [host.R_200,0.]
# initial_momentum = [-1.0*host.v_200,0.75*host.v_200]
initial_momentum = [0.8*host.v_200,0.9*host.v_200] # positive or negative - where should the minus sign be? (!!)
# initial_momentum = [0.,0.]

subhalo = Subhalo(mass_ratio=0.001,initial_position=initial_position,initial_momentum=initial_momentum)
# subhalo class might be unnecessary - combine with position/momentum classes? (!!)

print 'subhalo:'
print 'initial position = ', initial_position[0], initial_position[1]
print 'initial momentum = ', initial_momentum[0], initial_momentum[1]

position = ParticlePosition(subhalo.position)
momentum = ParticleMomentum(subhalo.momentum)

time = 0.
times = [time]
dt = 1.e14
while time < 1.e18:
	print time
	# print position.current_position
	force = gravity._calculate_gravitational_force(position=position.current_position)
	position._update_step(momentum=momentum.current_momentum,dt=dt)
	momentum._update_step(position=position.current_position,force=force,dt=dt)
	time+=dt
	times.append(time)
	# print time, position.position

times = np.array(times)
positions = np.array(position.position_array)
momenta = np.array(momentum.momentum_array)
print 'r min = ', positions[:,0].min()/host.R_200
print 'r max = ', positions[:,0].max()/host.R_200
print 'final time = ', time


forces = np.array(gravity.force_array)

plotting = 1
if plotting:
	plot = Plotting(times=times,positions=positions,momenta=momenta,forces=forces,host=host)
	plot.orbit()
	plot.orbit_color()
	plot.radial_position()
	plot.angular_position()
	plot.gravitational_force()
	plot.radial_velocity()



