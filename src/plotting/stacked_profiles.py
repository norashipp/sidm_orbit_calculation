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
# host = np.arange(51) # correct number?

integrator = 'leapfrog'
potential = 'spherical_NFW'
dt = 4e-3

v_thresh = 100 # km/s

subcount = 0

nbins = 13
rbins = np.linspace(0,3,nbins)
dr = rbins[1]-rbins[0]
rbc = rbins[1:]-dr/2

n = np.zeros_like(rbins[:-1])
nmt = np.zeros_like(rbins[:-1])
drag = 1
if drag:
	nd = np.zeros_like(rbins[:-1])

nsubs = np.ones_like(n)
nsubs_mt = np.ones_like(nmt)
nsubs_d = np.ones_like(nd)

sigma = np.ones_like(n)
sigma_mt = np.ones_like(nmt)
sigma_d = np.ones_like(nd)

var = np.ones_like(n)
var_mt = np.ones_like(nmt)
var_d = np.ones_like(nd)

for host_idx in hosts:
	infile = HOMEDIR+'sidm_orbit_calculation/src/output/final_positions_%i_%s_%s_%.0e.txt' %(host_idx,integrator,potential,dt)

	x,y,z = np.loadtxt(infile,unpack=True)
	r = np.sqrt(x**2 + y**2 + z**2)

	print 'Host %i' %host_idx
	print 'Number of subhalos = %i' %len(r)

	host = HostHalo(host_idx,potential)
	host.update(host.cosmo.age(0))

	subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)

	if drag:
		integrator = 'dissipative'
		infile = HOMEDIR+'sidm_orbit_calculation/src/output/final_positions_%i_%s_%s_%.0e.txt' %(host_idx,integrator,potential,dt)
		x,y,z = np.loadtxt(infile,unpack=True)
		rd = np.sqrt(x**2 + y**2 + z**2)
	
	ri = 0
	for i in range(len(host.subhalos)): # make this the outer loop, make list of subhalos to include, then bin them
		if host.subhalos[i]:
			# print i
			zz = 1/subs[i].a - 1
			tt = host.cosmo.age(zz)
			vmax = subs[i].v_max[-1]
			if vmax > v_thresh:
				subcount += 1
				# merger tree
				rf = np.sqrt(subs[i].rel_x[-1]**2 + subs[i].rel_y[-1]**2 + subs[i].rel_z[-1]**2)
				diff = np.abs(rf - rbc*host.R)
				rbin = diff.argmin()
				nmt[rbin] += 1
				# spherical NFW
				diff = np.abs(r[ri] - rbc*host.R)
				rbin = diff.argmin()
				n[rbin] += 1
				# drag force
				if drag:
					diff = np.abs(rd[ri] - rbc*host.R)
					rbin = diff.argmin()
					nd[rbin] += 1
				ri+=1

	vshell = 4/3*np.pi*((rbins[1:]*host.R)**3-(rbins[:-1]*host.R)**3)
	sd = nsubs*(dr*host.R)/vshell
	sd_mt = nsubs_mt*(dr*host.R)/vshell
	sd_d = nsubs_d*(dr*host.R)/vshell

	sigma*=sd/host.R
	sigma_mt*=sd_mt/host.R
	sigma_d*=sd_d/host.R

	var+=nsubs # /host.R? # error propagation!!
	var_mt+=nsubs_mt
	var_d+=nsubs_d

sigma = sigma**(1/len(hosts))
sigma_mt = sigma_mt**(1/len(hosts))
sigma_d = sigma_d**(1/len(hosts))

print sigma*host.R
print sigma_mt*host.R
print sigma_d*host.R

var/=len(hosts)
var_mt/=len(hosts)
var_d/=len(hosts)

err = np.sqrt(var)
err_mt = np.sqrt(var_mt)
err_d = np.sqrt(var_d)

# PLOTTING
plt.figure()
err = 0
if err:
	plt.errorbar(rbc,sigma,xerr=dr,yerr=err,ls='-',color='b',label=r'$\mathrm{Spherical\ NFW}$',lw=3,markersize='10')
	plt.errorbar(rbc,sigma_d,xerr=dr,yerr=err_d,ls='-',color='r',label=r'$\mathrm{Drag\ Force}$',lw=3,markersize='10')
	plt.errorbar(rbc,sigma_mt,xerr=dr,yerr=err_mt,ls='-',color='g',label=r'$\mathrm{Merger\ Tree}$',lw=3,markersize='10')
else:
	plt.plot(rbc,sigma,'-',color='b',label=r'$\mathrm{Spherical\ NFW}$',lw=3,markersize='10')
	plt.plot(rbc,sigma_d,'-',color='r',label=r'$\mathrm{Drag\ Force}$',lw=3,markersize='10')
	plt.plot(rbc,sigma_mt,'-',color='g',label=r'$\mathrm{Merger\ Tree}$',lw=3,markersize='10')
plt.grid()
plt.xlabel(r'$\mathrm{r/R_{200m}}$')
plt.ylabel(r'$\mathrm{\Sigma /R_{200m}\ (subhalos/Mpc^3)}$')
plt.title(r'$\mathrm{Stacked,\ v_{thresh}\ =\ %.2f\ km/s}$' %(v_thresh))
plt.legend()
plt.yscale('log')
# plt.savefig('plots/subhalo_distribution_%i_%s_%s_%.0e_%i.png'  %(host_idx,integrator,potential,dt,v_thresh))
plt.show()
