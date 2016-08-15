#!/usr/bin/env python

from __future__ import division
import numpy as np
import cPickle
import os
import sys
import cProfile
from time import sleep

from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.halos.subhalo import *
from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.calculation.integrate import *
from sidm_orbit_calculation.src.utils.setup import *

class Sim:
	# def __init__(self, host_halo_mass, host_idx, subhalo_mass, dt, tmax, integration_method, potential, initial_position, intiial_momentum):
	def __init__(self, host_idx, dt, integration_method, potential):
		self.dt = dt
		# self.tmax = tmax

		# self.host = HostHalo(M=host_halo_mass, potential=potential, idx=host_idx, tmax=tmax)
		self.host = HostHalo(idx=host_idx, potential=potential)
		
		self.integrate = integrator_dict[integration_method]

		# self.initiate_subhalo(mass=subhalo_mass, initial_position=initial_position, initial_momentum=initial_momentum)

	# def initiate_subhalo(self, mass, initial_position, initial_momentum):
	# 	self.subhalo = Subhalo(host=self.host, M=mass, initial_position=initial_position, initial_momentum=initial_momentum)

	def sim(self, subhalo, printing=0, writing=0, outfile='pickle.dat'):
		self.time = subhalo.t0
		self.host.update(self.time)
		times = [self.time]

		while self.time < self.host.cosmo.age(0): # assuming all subhalos survive to z=0
			self.host.update(self.time)
			self.integrate(subhalo=subhalo,dt=self.dt)
			self.time+=self.dt
			
			# for now, just assume always using merger tree output - can add in class inheritance later
			# if self.host.host_idx != None:
			# self.host.update(self.time)

			subhalo.calculate_energy() # weird - reorganize

			# if not subhalo.count % 500:
			times.append(np.copy(self.time))
			
			subhalo.count += 1

		times = np.array(times)
		positions = np.asarray(subhalo.position_array)
		momenta = np.asarray(subhalo.momentum_array)
		gravity = np.asarray(subhalo.gravity.force_array)
		drag = np.asarray(subhalo.drag.force_array)
		density = np.asarray(self.host.density_array)
		energy = np.asarray(subhalo.energy_array)

		# self.host.M_sp = None
		# self.host.R_s_sp = None
		# self.host.q_sp = None
		# self.host.s_sp = None
		# self.host.mass_function = None
		# self.host.density_function = None
		# self.host.potential_function = None

		'''
		d = dir(self.host)
		for x in d: print x, eval("type(self.host.%s)" % x)
		print

		d = dir(self.host.subhalos[0])
		for x in d: print x, eval("type(self.host.subhalos[0].%s)" % x)
		'''

		self.output = [times,positions,momenta,gravity,drag,density,energy,self.host.host_idx,self.host.potential,host.R]

		if writing:
			self.write_output(outfile)

		if printing:
			print 'r min = ', positions[:,0].min()/self.host.R
			print 'r max = ', positions[:,0].max()/self.host.R
			print 'r initial = ', positions[:,0][0]/self.host.R
			print 'r final = ', positions[:,0][-1]/self.host.R
			print 'final time = ', times.max()
			print 'final gravitational force = %10.5g %10.5g' % (subhalo.gravity.vector(phi)[0],subhalo.gravity.vector(phi)[1])
			print 'final drag force = %10.5g %10.5g' % (subhalo.drag.force[0],subhalo.drag.force[1])

	def write_output(self,output_file):
		f = open(output_file,'wb')
		cPickle.dump(self.output,f)
		f.close()

###########################################

host_idx = int(sys.argv[1])
dt = float(sys.argv[2])# /1e9 # input in Gyr now!
integrator = sys.argv[3]
potential = sys.argv[4]

'''
host_idx = int(sys.argv[1])
host_halo_mass = float(sys.argv[2])
subhalo_mass = float(sys.argv[3])
dt = float(sys.argv[4])# /1e9 # input in Gyr now!
tmax = float(sys.argv[5])# /1e9
integrator = sys.argv[6]
potential = sys.argv[7]
index = int(sys.argv[8])
# print len(sys.argv)

try:
    initial_position = np.array([float(sys.argv[9]), float(sys.argv[10]), float(sys.argv[11])])
    initial_momentum = np.array([float(sys.argv[12]), float(sys.argv[13]), float(sys.argv[14])])
except:
    initial_position = np.array([0,0,0])
    initial_momentum = np.array([0,0,0])
'''

my_sim = Sim(host_idx, dt, integrator, potential)
# my_sim = Sim(1e14, 1e12, 1e4/seconds_to_years, 1e10/seconds_to_years, 'leapfrog', 'point_mass', np.array([1e22,0,0]), np.array([0,0,0]))

t0 = my_sim.host.cosmo.age(0)
for i in range(len(my_sim.host.subhalos)):
	subhalo = my_sim.host.subhalos[i]
	sub_idx = i
	print 'Integrating subhalo %i/%i' %(sub_idx, len(my_sim.host.subhalos))
	outfile = HOMEDIR + 'sidm_orbit_calculation/src/output/%i_%s_%s_%.2e_%i.dat' %(host_idx, integrator, potential, dt, sub_idx)
	my_sim.sim(subhalo=subhalo, printing=0, writing=1, outfile=outfile)

# example: python sim.py 1e14 1e12 1e4 1e10 leapfrog spherical_NFW 0
# cProfile.run('my_sim.sim(printing=0, writing=1, plotting=0, outfile=outfile)')
