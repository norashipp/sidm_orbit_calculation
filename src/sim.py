from __future__ import division
import numpy as np
import cPickle

from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
# from sidm_orbit_calculation.src.timestep.particle_evolution import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.halos.subhalo import *
from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.calculation.integrate import *

class Sim:
	def __init__(self,dt=1.e12,tmax=1.e14,integration_method='leapfrog',potential='point_mass'):
		self.host = HostHalo(M=1.e13*M_sol,potential=potential)
		self.dt = dt
		self.tmax = tmax
		self.time = 0.
		self.orbit_count = 0

		self.integrate = integrator_dict[integration_method]

		initial_position, initial_momentum = self.initial_parameters()
		self.initiate_subhalo(position=initial_position,momentum=initial_momentum)
		
	def initial_parameters(self):
		# FOR NOW MANUALLY DEFINE INITIAL PARAMETERS
		initial_position = np.array([self.host.R_200,0,0])
		initial_momentum = np.array([0.,3e5,0.])
		return initial_position, initial_momentum

	def initiate_subhalo(self,position,momentum):
		self.subhalo = Subhalo(host=self.host,mass_ratio=0.5,initial_position=position,initial_momentum=momentum)
		
	def sim(self,printing=0,writing=0,plotting=0):
		times = [self.time]
		
		while self.time < self.tmax:
			self.integrate(subhalo=self.subhalo,dt=self.dt)
			self.time+=self.dt
			times.append(self.time)
			# phi = np.arctan2(self.subhalo.position[1],self.subhalo.position[0])
			# if printing == 2:
			# 	print 'time = ', self.time
			# 	print 'position = %10.5g %10.5g' % (self.subhalo.position[0],self.subhalo.position[1])
			# 	print 'force = %10.5g %10.5g' % (self.subhalo.gravity.vector(phi)[0],self.subhalo.gravity.vector(phi)[1])
				# assert 0
			# self.orbit_count+=1

		times = np.array(times)
		positions = np.asarray(self.subhalo.position_array)
		momenta = np.asarray(self.subhalo.momentum_array)
		gravity = np.asarray(self.subhalo.gravity.force_array)
		drag = np.asarray(self.subhalo.drag.force_array)
		density = np.asarray(self.host.density_array)
		self.output = [times,positions,momenta,gravity,drag,density,self.host]

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
			print 'final gravitational force = %10.5g %10.5g' % (self.subhalo.gravity.vector(phi)[0],self.subhalo.gravity.vector(phi)[1])
			print 'final drag force = %10.5g %10.5g' % (self.subhalo.drag.force[0],self.subhalo.drag.force[1])

	def write_output(self):
		f = open('data/pickle.dat','wb')
		cPickle.dump(self.output,f)
		f.close()

	def plot_output(self):
		times,positions,momenta,gravity,drag,density,host = self.output
		plot = Plotting(times=times,positions=positions,momenta=momenta,gravity=gravity,drag=drag,density=density,host= host)
		plot.orbit()
		# plot.drag_force()
		# plot.density()
		# plot.test_density()
		# plot.orbit_color()
		# plot.radial_position()
		# plot.angular_position()
		plot.gravitational_force()
		# plot.radial_velocity()

dt = 1e4/seconds_to_years
#simulation = Sim(dt=dt,tmax=2e4*dt,integration_method='euler')
simulation = Sim(dt=dt,tmax=1e10/seconds_to_years,integration_method='leapfrog',potential='point_mass') # dynamical time 5e9 years
simulation.sim(printing=0,writing=1,plotting=0) # for max printing = 2

