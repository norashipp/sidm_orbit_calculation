from __future__ import division
import numpy as np
import cPickle

from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.timestep.particle_evolution import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.halos.subhalo import *
from sidm_orbit_calculation.src.plotting.make_plots import *

class Sim:
	def __init__(self,dt=1.e12,tmax=1.e14,integration_method='leapfrog',potential='point_mass'):
		self.host = HostHalo(M=1.e13*M_sol,potential=potential)
		self.gravity = GetGravitationalForce(self.host)
		self.dt = dt
		self.tmax = tmax
		self.time = 0.
		self.loop_count = 0

		self.integration_method = integration_method

		initial_position, initial_momentum = self.initial_parameters()
		self.initiate_particle(position=initial_position,momentum=initial_momentum,integration_method=self.integration_method)
		
	def initial_parameters(self):
		# FOR NOW MANUALLY DEFINE INITIAL PARAMETERS
		r0 = self.host.R_200
		phi0 = 0.
		initial_position = np.array([r0*np.cos(phi0),r0*np.sin(phi0)]) # cartesian coordinates
		initial_momentum = np.array([0.,2.2e5])
		# initial_momentum = np.array([0.,0.])
		return initial_position, initial_momentum

	def initiate_particle(self,position=None,momentum=None,integration_method='leapfrog'):
		self.position = ParticlePosition(initial_position=position,dt=dt)
		self.momentum = ParticleMomentum(integration_method=integration_method,initial_momentum=momentum,gravity=GetGravitationalForce(self.host),dt=dt,position=position)
		
	def step(self):
		# probably get rid of this if statement at some point... (!!)
		if self.integration_method in ['euler', 'leapfrog']:
			if self.integration_method == 'euler':
				current_position = self.position.current_position[:]
			self.position.update_step(momentum=self.momentum.current_momentum)
			if self.integration_method == 'leapfrog':
				current_position = self.position.current_position[:]
			self.momentum.update_step(position=current_position)
		
	def sim(self,printing=0,writing=0,plotting=0):
		times = [self.time]
		phi = 0.
		while self.time < self.tmax:
			self.step()
			self.time+=self.dt
			times.append(self.time)
			phi = np.arctan(self.position.current_position[1]/self.position.current_position[0])
			if printing == 2:
				print 'time = ', self.time
				print 'position = %10.5g %10.5g' % (self.position.current_position[0],self.position.current_position[1])
				phi = np.arctan(self.position.current_position[1]/self.position.current_position[0])
				print 'force = %10.5g %10.5g' % self.momentum.gravity.vector(phi)
				# assert 0
		times = np.array(times)
		positions = np.array(self.position.position_array)
		momenta = np.array(self.momentum.momentum_array)
		forces = np.array(self.momentum.gravity.force_array)
		self.output = [times,positions,momenta,forces]

		if writing:
			self.write_output()

		if plotting:
			self.plot_output()

		if printing:
			print 'r min = ', positions[:,0].min()/self.host.R_200
			print 'r max = ', positions[:,0].max()/self.host.R_200
			print 'r initial = ', positions[:,0][0]/self.host.R_200
			print 'r final = ', positions[:,0][-1]/self.host.R_200
			print 'final time = ', times.max()
			print 'final force = %10.5g %10.5g' % self.momentum.gravity.vector(phi)

	def write_output(self):
		f = open('data/pickle.dat','wb')
		cPickle.dump(self.output,f)
		f.close()

	def plot_output(self):
		times,positions,momenta,forces = self.output
		plot = Plotting(times=times,positions=positions,momenta=momenta,forces=forces,host=self.host)
		plot.orbit()
		# plot.orbit_color()
		# plot.radial_position()
		# plot.angular_position()
		# plot.gravitational_force()
		# plot.radial_velocity()

dt = 4e7/seconds_to_years
#simulation = Sim(dt=dt,tmax=2e4*dt,integration_method='euler')
simulation = Sim(dt=dt,tmax=10 * 5e9/seconds_to_years,integration_method='leapfrog',potential='spherical_NFW')
simulation.sim(printing=1,writing=1,plotting=1) # for max printing = 2

