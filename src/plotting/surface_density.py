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

v_thresh = 0 # km/s

for host_idx in hosts:
	integrator = 'leapfrog'
	infile = HOMEDIR+'sidm_orbit_calculation/src/output/final_positions/final_positions_%i_%s_%s_%.0e.txt' %(host_idx,integrator,potential,dt)

	pos = np.loadtxt(infile)
	r = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
	print 'Number of subhalos = %i' %len(r)

	host = HostHalo(host_idx,potential)
	host.update(host.cosmo.age(0))

	subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)

	subcount = 0

	nbins = 13
	rbins = np.linspace(0,3*host.R,nbins)
	dr = rbins[1]-rbins[0]
	rbc = rbins[1:]-dr/2
	# print rbins
	nsubs = np.zeros_like(rbins[:-1])
	nsubs_mt = np.zeros_like(rbins[:-1])

	drag = 1
	if drag:
		integrator = 'dissipative'
		infile = HOMEDIR+'sidm_orbit_calculation/src/output/final_positions/final_positions_%i_%s_%s_%.0e.txt' %(host_idx,integrator,potential,dt)
		pos = np.loadtxt(infile)
		rd = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
		nsubs_d = np.zeros_like(rbins[:-1])
	
	vmax_0 = []
	vmax_max = []

	ri = 0
	for i in range(len(host.subhalos)): # make this the outer loop, make list of subhalos to include, then bin them
		if host.subhalos[i]:
			# print i
			zz = 1/subs[i].a - 1
			tt = host.cosmo.age(zz)
			# vmax_sp = interp1d(tt,subs[i].v_max)
			# vmax_acc = vmax_sp(host.subhalos[i].t0)
			
			# vmax_vals = subs[i].v_max
			# plt.plot(tt,vmax_vals)

			vmax = subs[i].v_max[-1]

			if vmax > v_thresh:
				subcount += 1

				# determine radial bins for each

				# merger tree
				rf = np.sqrt(subs[i].rel_x[-1]**2 + subs[i].rel_y[-1]**2 + subs[i].rel_z[-1]**2)
				'''
				diff = rf - rbins[:-1]
				diff[diff<0] = np.inf
				rbin = diff.argmin()
				'''
				diff = np.abs(rf-rbc)
				rbin = diff.argmin()
				nsubs_mt[rbin] += 1
				
				# spherical NFW
				'''
				diff = r[ri] - rbins[:-1]
				diff[diff<0] = np.inf
				'''
				diff = np.abs(r[ri] - rbc)
				rbin = diff.argmin()
				nsubs[rbin] += 1

				# drag force
				if drag:
					'''
					diff = rd[ri] - rbins[:-1]
					diff[diff<0] = np.inf
					'''
					diff = np.abs(rd[ri] - rbc)
					rbin = diff.argmin()
					nsubs_d[rbin] += 1

				ri+=1
			# else:
			# 	print 'not updating'
        # else:
        # 	print 'skipping %i' %i

	vshell = 4/3*np.pi*(rbins[1:]**3-rbins[:-1]**3)
	sd = nsubs*dr/vshell
	sd_mt = nsubs_mt*dr/vshell
	sd_d = nsubs_d*dr/vshell

	print np.sum(nsubs)

	# print nsubs, np.sqrt(nsubs)
	# print nsubs_d, np.sqrt(nsubs_d)
	# print nsubs_mt, np.sqrt(nsubs_mt)

	# print 'results'
	# print nsubs, np.sum(nsubs)
	# print nsubs_mt, np.sum(nsubs_mt)
	# print subcount

	rbc/=host.R

	# PLOTTING
	plt.figure()
	# plt.bar(rbins[:-1],nsubs/vshell,width=rbins[1]-rbins[0],color='c',alpha=0.8,label=r'$\mathrm{Spherical\ NFW}$',zorder=1)
	# plt.bar(rbins[:-1],nsubs_mt/vshell,width=rbins[1]-rbins[0],color='g',alpha=1.0,label=r'$\mathrm{Merger\ Tree}$',zorder=0)
	# plt.bar(rbins[:-1],nsubs_d/vshell,width=rbins[1]-rbins[0],color='m',alpha=0.5,label=r'$\mathrm{Drag\ Force}$',zorder=2)
	err = 0
	if err:
		plt.errorbar(rbc,sd,xerr=dr,yerr=np.sqrt(nsubs),ls='-',color='b',label=r'$\mathrm{Spherical\ NFW}$',lw=3,markersize='10')
		plt.errorbar(rbc,sd_d,xerr=dr,yerr=np.sqrt(nsubs_d),ls='-',color='r',label=r'$\mathrm{Drag\ Force}$',lw=3,markersize='10')
		plt.errorbar(rbc,sd_mt,xerr=dr,yerr=np.sqrt(nsubs_mt),ls='-',color='g',label=r'$\mathrm{Merger\ Tree}$',lw=3,markersize='10')
	else:
		plt.plot(rbc,sd,'-',color='b',label=r'$\mathrm{Spherical\ NFW}$',lw=3,markersize='10')
		plt.plot(rbc,sd_d,'--',color='r',label=r'$\mathrm{Drag\ Force}$',lw=3,markersize='10')
		plt.plot(rbc,sd_mt,'-',color='g',label=r'$\mathrm{Merger\ Tree}$',lw=3,markersize='10')
	plt.grid()
	plt.xlabel(r'$\mathrm{r/R_{200m}}$')
	plt.ylabel(r'$\mathrm{Surface\ Density\ (subhalos/Mpc^2)}$')
	plt.title(r'$\mathrm{Host\ %i,\ v_{thresh}\ =\ %.2f\ km/s}$' %(host_idx,v_thresh))
	plt.legend()
	plt.yscale('log')
	plt.xscale('log')
	plt.xlim(0,3*host.R)
	# plt.savefig('plots/subhalo_distribution_%i_%s_%s_%.0e_%i.png'  %(host_idx,integrator,potential,dt,v_thresh))
plt.show()
