#!/usr/bin/env python

from __future__ import division
import numpy as np
import cPickle
import os
import sys

from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
# from sidm_orbit_calculation.src.timestep.particle_evolution import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.halos.subhalo import *
from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.calculation.integrate import *

class Sim:
	def __init__(self, host_halo_mass, subhalo_mass, dt, tmax, integration_method, potential):
		self.host = HostHalo(M=host_halo_mass*M_sol, potential=potential)
		self.dt = dt
		self.tmax = tmax
		self.time = 0.
		self.orbit_count = 0

		self.integrate = integrator_dict[integration_method]

		# initial_position, initial_momentum = self.initial_parameters()
		self.initiate_subhalo(mass=subhalo_mass,position=None,momentum=None)
		
	# def initial_parameters(self):
	# 	initial_position, initial_momentum = initial_conditions(self.subhalo)
	# 	initial_position = np.array([self.host.R_200,0,0])
	# 	initial_momentum = np.array([0.,3e5,0.])
	# 	return initial_position, initial_momentum

	def initiate_subhalo(self, mass, position, momentum):
		self.subhalo = Subhalo(host=self.host, M=mass, initial_position=position, initial_momentum=momentum)
		
	def sim(self, printing=0, writing=0, plotting=0, outfile='pickle.dat'):
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
			self.write_output(outfile)

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

	def write_output(self,output_file):
		# fn = '/Users/nora/sidm_orbit_calculation/src/data/' + str(self.host.M) + '_' + str(self.subhalo.M) + '_'
		# i = 0
		# while os.path.isfile(fn+str(i)+'.dat'):
		# 	i+=0
		# fname = fn + str(i) + '.dat'
		
                f = open(output_file,'wb')
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

###########################################

# dt = 1e4/seconds_to_years
#simulation = Sim(dt=dt,tmax=2e4*dt,integration_method='euler')
# simulation = Sim(dt=dt,tmax=2e10/seconds_to_years,integration_method='leapfrog',potential='spherical_NFW') # dynamical time 5e9 years
# simulation.sim(printing=0,writing=1,plotting=0) # for max printing = 2

###########################################

host_halo_mass = float(sys.argv[1])
subhalo_mass = float(sys.argv[2])
dt = float(sys.argv[3])/seconds_to_years
tmax = float(sys.argv[4])/seconds_to_years
integrator = sys.argv[5]
potential = sys.argv[6]
index = int(sys.argv[7])

homedir = '/home/norashipp/sidm_orbit_calculation/'
outfile = homedir + 'src/output/%.1e_%.1e_%.1e_%.1e_%s_%s_%i.dat' %(host_halo_mass, subhalo_mass, dt*seconds_to_years, tmax*seconds_to_years, integrator, potential, index)

my_sim = Sim(host_halo_mass, subhalo_mass, dt, tmax, integrator, potential)
my_sim.sim(printing=0, writing=1, plotting=0, outfile=outfile)

# maybe make it easier to not include all of these parameters and default to typical values?

# example: python sim.py 1e14 1e12 1e4 1e10 leapfrog spherical_NFW 0
