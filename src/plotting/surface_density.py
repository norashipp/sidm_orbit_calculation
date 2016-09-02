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

dpi = 175
fontsize = 15
plt.rc('savefig', dpi=dpi)
plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.rc('xtick.major', pad=5)
plt.rc('xtick.minor', pad=5)
plt.rc('ytick.major', pad=5)
plt.rc('ytick.minor', pad=5)

hosts = np.array(sys.argv[1:],dtype=int)
# host_idx = int(sys.argv[1])
integrator = 'leapfrog'
potential = 'spherical_NFW'
dt = 4e-3

v_thresh = 150 # km/s

for host_idx in hosts:
	infile = HOMEDIR+'sidm_orbit_calculation/src/output/final_positions_%i_%s_%s_%.0e.txt' %(host_idx,integrator,potential,dt)

	pos = np.loadtxt(infile)
	r = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
	print 'Number of subhalos = %i' %len(r)

	host = HostHalo(host_idx,potential)
	host.update(host.cosmo.age(0))

	subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)

	subcount = 0

	nbins = 13
	rbins = np.linspace(0,3*host.R,nbins)
	# print rbins
	nsubs = np.zeros_like(rbins[:-1])
	nsubs_mt = np.zeros_like(rbins[:-1])
	ri = 0
	for i in range(len(host.subhalos)): # make this the outer loop, make list of subhalos to include, then bin them
		if host.subhalos[i]:
			# print i
			zz = 1/subs[i].a - 1
			tt = host.cosmo.age(zz)
			vmax_sp = interp1d(tt,subs[i].v_max)
			vmax_acc = vmax_sp(host.subhalos[i].t0)
			if vmax_acc > v_thresh:
				subcount += 1
				# determine radial bins for each
				# merger tree
				rf = np.sqrt(subs[i].rel_x[-1]**2 + subs[i].rel_y[-1]**2 + subs[i].rel_z[-1]**2)
				diff = rf - rbins[:-1]
				diff[diff<0] = np.inf
				rbin = diff.argmin()
				# print 'merger tree'
				# print rf
				# print rbin
				nsubs_mt[rbin] += 1

				# calculation
				diff = r[ri] - rbins[:-1]
				ri+=1
				diff[diff<0] = np.inf
				rbin = diff.argmin()
				# print 'calculation'
				# print r[i]
				# print rbin
				# print 
				nsubs[rbin] += 1
		# else:
			# print 'skipping %i' %i
	
	vshell = 4/3*np.pi*(rbins[1:]**3-rbins[:-1]**3)

	# print 'results'
	# print nsubs, np.sum(nsubs)
	# print nsubs_mt, np.sum(nsubs_mt)
	# print subcount

	plt.figure()
	plt.bar(rbins[:-1],nsubs/vshell,width=rbins[1]-rbins[0],color='c',alpha=0.5,label=r'$\mathrm{Spherical\ NFW}$')
	plt.bar(rbins[:-1],nsubs_mt/vshell,width=rbins[1]-rbins[0],color='g',alpha=0.5,label=r'$\mathrm{Merger\ Tree}$',zorder=0)
	# plt.hist(np.ones_like(rbins[:-1]), weights=nsubs, color='c', bins=rbins, histtype='step', lw=3, normed=True)
	# plt.hist(np.ones_like(rbins[:-1]), weights=nsubs_mt, color='g', bins=rbins, histtype='step', lw=3, normed=True)
	plt.grid()
	plt.xlabel(r'$\mathrm{x/R_{200m}}$')
	plt.ylabel(r'$\mathrm{Surface\ Density\ (subhalos/Mpc^3)}$')
	plt.title(r'$\mathrm{Host\ %i,\ v_{thresh}\ =\ %.2f\ km/s}$' %(host_idx,v_thresh))
	plt.legend()
	plt.yscale('log')
	plt.savefig('plots/subhalo_distribution_%i_%s_%s_%.0e_%i.png'  %(host_idx,integrator,potential,dt,v_thresh))

# plt.show()
