from __future__ import division
import numpy as np
import cPickle

from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.timestep.particle_evolution import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.halos.subhalo import *
from sidm_orbit_calculation.src.plotting.make_plots import *

host = HostHalo(M=1.e13*M_sol,potential='point_mass')
gravity = GetGravitationalForce(host) # will append to same array for all particles, need to change this if plottings

def step(time=None,dt=1.e12):
	force_r,force_phi,dxdiv = gravity._calculate_gravitational_force(position=position.current_position)
	force = (force_r,force_phi)
	position.update_step(momentum=momentum.current_momentum,dt=dt)
	momentum.update_step(position=position.current_position,force=force,dt=dt)
	return None

def initiate_particles(n_particles):
	particle_positions = []
	particle_momenta = []
	for i in range(n_particles):
		paritcle_positions.append(ParticlePosition(initial_position))
		particle_momenta.append(ParticleMomentum(initial_momentum))
	return particle_positions, particle_momenta

def sim(position=None,momentum=None,dt=1.e12):
	time = 0.
	times = [time]
	while position.current_position[1] < 20*(2*np.pi):
		step(time,dt)
		time+=dt
		times.append(time)
	times = np.array(times)
	positions = np.array(position.position_array)
	momenta = np.array(momentum.momentum_array)
	forces = np.array(gravity.force_array)
	return times,positions,momenta,forces

def save_output(data=None,indx=None):
	fname = 'pickle'+str(indx)+'.dat'
	ff = file.open(fname)
	cPickle.dump(data)
	ff.close()

n_particles = 10
initial_position = [host.R_200,0.]
initial_momentum = [-0.8*host.v_200,0.9*host.v_200]
particle_positions,particle_momenta = initiate_particles(initial_position,initial_momentum)

for ii in range(n_paritlces):
	times,positions,momenta,forces = sim(position=paritcle_positions[ii],momentum=particle_momenta[ii],dt = 1.e12)
	save_output(data=[times,positions,momenta,forces],indx=ii)




plotting = 1
if plotting:
	plot = Plotting(times=times,positions=positions,momenta=momenta,forces=forces,host=host)
	plot.multiple_orbits()
	