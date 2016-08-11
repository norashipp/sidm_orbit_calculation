#!/usr/bin/env python

from __future__ import division
import numpy as np
import cPickle
import os
import sys
import cProfile

from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.halos.subhalo import *
from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.calculation.integrate import *
from sidm_orbit_calculation.src.utils.setup import *
import sidm_orbit_calculation.src.merger_tree.cluster as cluster

class Sim:
	def __init__(self, host_halo_mass, subhalo_mass, dt, tmax, integration_method, potential, initial_position, intiial_momentum):
		
		self.initiate_host_halo(M=host_halo_mass, potential=potential)

		self.dt = dt
		self.tmax = tmax
		self.time = 0.
		self.orbit_count = 0

		self.integrate = integrator_dict[integration_method]

		self.initiate_subhalo(mass=subhalo_mass, initial_position=initial_position, initial_momentum=initial_momentum)
	
	def initiate_host_halo(self, M, potential):
		if M:
			self.host = HostHalo(M=M, potential=potential)
		else:
			print 'Importing host halo parameters from merger tree...'
			hs = cluster.HostHalos(homedir + 'sidm_orbit_calculation/src/merger_tree/clusters.dat')
			host_idx = 40
			snap = 0

			a = hs[host_idx].a[snap] 
			M = hs[host_idx].m_200m[snap]
			R_s = hs[host_idx].r_s[snap]
			q = hs[host_idx].b_to_a[snap]
			s = hs[host_idx].c_to_a[snap]

			self.host = HostHalo(M=M, potential=potential, s = s, q = q, a = a, R_s = R_s)

	def initiate_subhalo(self, mass, initial_position, initial_momentum):
		self.subhalo = Subhalo(host=self.host, M=mass, initial_position=initial_position, initial_momentum=initial_momentum)

	def sim(self, printing=0, writing=0, outfile='pickle.dat'):
		times = [self.time]
		
		while self.time < self.tmax:
			self.integrate(subhalo=self.subhalo,dt=self.dt)
			self.time+=self.dt

			self.subhalo.calculate_energy()

			# if not self.subhalo.count % 500:
			times.append(self.time)
			self.subhalo.count += 1

		times = np.array(times)
		positions = np.asarray(self.subhalo.position_array)
		momenta = np.asarray(self.subhalo.momentum_array)
		gravity = np.asarray(self.subhalo.gravity.force_array)
		drag = np.asarray(self.subhalo.drag.force_array)
		density = np.asarray(self.host.density_array)
		energy = np.asarray(self.subhalo.energy_array)
		self.output = [times,positions,momenta,gravity,drag,density,energy,self.host]

		if writing:
			self.write_output(outfile)

		if printing:
			print 'r min = ', positions[:,0].min()/self.host.R
			print 'r max = ', positions[:,0].max()/self.host.R
			print 'r initial = ', positions[:,0][0]/self.host.R
			print 'r final = ', positions[:,0][-1]/self.host.R
			print 'final time = ', times.max()
			print 'final gravitational force = %10.5g %10.5g' % (self.subhalo.gravity.vector(phi)[0],self.subhalo.gravity.vector(phi)[1])
			print 'final drag force = %10.5g %10.5g' % (self.subhalo.drag.force[0],self.subhalo.drag.force[1])

	def write_output(self,output_file):
		f = open(output_file,'wb')
		cPickle.dump(self.output,f)
		f.close()

###########################################

host_halo_mass = float(sys.argv[1])
subhalo_mass = float(sys.argv[2])
dt = float(sys.argv[3])# /1e9 # input in Gyr now!
tmax = float(sys.argv[4])# /1e9
integrator = sys.argv[5]
potential = sys.argv[6]
index = int(sys.argv[7])
# print len(sys.argv)

try:
    initial_position = np.array([float(sys.argv[8]), float(sys.argv[9]), float(sys.argv[10])])
    initial_momentum = np.array([float(sys.argv[11]), float(sys.argv[12]), float(sys.argv[13])])
except:
    initial_position = np.array([0,0,0])
    initial_momentum = np.array([0,0,0])

homedir = home_directory()
outfile = homedir + 'sidm_orbit_calculation/src/output/%.1e_%.1e_%.1e_%.1e_%s_%s_%i.dat' %(host_halo_mass, subhalo_mass, dt*seconds_to_years, tmax*seconds_to_years, integrator, potential, index)

my_sim = Sim(host_halo_mass, subhalo_mass, dt, tmax, integrator, potential, initial_position, initial_momentum)
# my_sim = Sim(1e14, 1e12, 1e4/seconds_to_years, 1e10/seconds_to_years, 'leapfrog', 'point_mass', np.array([1e22,0,0]), np.array([0,0,0]))

# cProfile.run('my_sim.sim(printing=0, writing=1, plotting=0, outfile=outfile)')
my_sim.sim(printing=0, writing=1, outfile=outfile)

# example: python sim.py 1e14 1e12 1e4 1e10 leapfrog spherical_NFW 0
