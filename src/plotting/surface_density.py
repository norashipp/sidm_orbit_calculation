from __future__ import division
import pylab
import numpy as np
import cPickle
import sys
import glob
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt

from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.utils.setup import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.merger_tree.cluster import *

hosts = np.array(sys.argv[1:],dtype=int)
# host_idx = int(sys.argv[1])
integrator = 'leapfrog'
potential = 'spherical_NFW'
dt = 4e-3

v_thresh = 150 # km/s

for host_idx in hosts:
	infile = HOMEDIR+'sidm_orbit_calculation/src/output/final_positions_%i_%s_%s_%.0e.txt' %(host_idx,integrator,potential,dt)

	pos = np.loadtxt(infile)
	dist = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)

	host = HostHalo(host_idx,potential)
	host.update(host.cosmo.age(0))

	subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)
	mass = []
	dist_mt = []

	nbins = 15
	nsubs = np.zeros_like(nbins)
	for rbin in range(nbins):
		for i in range(len(host.subhalos)): # make this the outer loop, make list of subhalos to include, then bin them
			if host.subhalos[i]:
				zz = 1/subs[i].a - 1
				tt = host.cosmo.age(zz)

				vmax_sp = interp1d(tt,subs[i].v_max)
				vmax_acc = v_max_sp(host.subhalos[i].t0)
				if vmax_acc > v_thresh:
					print vmax_acc
					nsubs[rbin] += 1
				
				# print i
				mass.append(subs[i].m_200m[-1])
				dist_mt.append(np.sqrt(subs[i].rel_x[-1]**2 + subs[i].rel_y[-1]**2 + subs[i].rel_z[-1]**2))
			# else:
				# print 'skipping %i' %i
		mass = np.asarray(mass)
		dist_mt = np.asarray(dist_mt)

	print host_idx
	print dist.shape, mass.shape, dist_mt.shape
	print 

	bins = np.linspace(0,5,21)
	plt.figure()
	plt.hist(dist/host.R, color='c', bins=bins, weights=mass, histtype='step', lw=3, normed=True)
	plt.hist(dist_mt/host.R, color='g', bins=bins, weights=mass, histtype='step', lw=3, normed=True)
	plt.grid()
	plt.xlabel('x/R200m')
	plt.ylabel('mass weighted subhalo distribution')
	plt.show()
	# plt.savefig('subhalo_distribution_%i_%s_%s_%.0e.png'  %(host_idx,integrator,potential,dt))
