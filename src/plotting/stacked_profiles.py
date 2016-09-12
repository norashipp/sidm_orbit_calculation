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
if len(hosts) == 0: hosts = np.arange(1,51)
nhosts = len(hosts)

integrator = 'leapfrog'
potential = 'spherical_NFW'
dt = 4e-3
drag = 1
triaxial = 1

v_thresh = 100 # km/s

nbins = 12
rbins = np.linspace(0,3,nbins+1)
dr = rbins[1]-rbins[0]
rbc = rbins[1:]-dr/2

# sigma = np.ones_like(nbins)
# sigma_mt = np.ones_like(nbins)
# sigma_d = np.ones_like(nbins)

sigma = np.ones((nhosts,nbins))
sigma_mt = np.ones((nhosts,nbins))
sigma_d = np.ones((nhosts,nbins))
sigma_t = np.ones((nhosts,nbins))

# var = np.zeros_like(nbins)
# var_mt = np.zeros_like(nbins)
# var_d = np.zeros_like(nbins)

for j in range(nhosts):
	host_idx = hosts[j]
	# spherical NFW
	infile = HOMEDIR+'Dropbox/SIDMdata/final_positions/final_positions_%i_leapfrog_spherical_NFW_%.0e.txt' %(host_idx,dt)
	x,y,z = np.loadtxt(infile,unpack=True)
	r = np.sqrt(x**2 + y**2 + z**2)

	host = HostHalo(host_idx,potential,subs=False)
	host.update(host.cosmo.age(0))
	print 'Host %i' %host_idx
	print 'M = %.2e' %host.M
	print 'Number of subhalos = %i' %len(r)
	
	subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)

	if drag:
		sig = 21
		infile = HOMEDIR+'Dropbox/SIDMdata/final_positions/sigma%i/final_positions_%i_dissipative_%s_%.0e.txt' %(sig,host_idx,potential,dt)
		x,y,z = np.loadtxt(infile,unpack=True)
		rd = np.sqrt(x**2 + y**2 + z**2)
		
	if triaxial:
		infile = HOMEDIR+'Dropbox/SIDMdata/final_positions/final_positions_%i_leapfrog_triaxial_NFW_BT_%.0e.txt' %(host_idx,dt)
		x,y,z = np.loadtxt(infile,unpack=True)
		rt = np.sqrt(x**2 + y**2 + z**2)

	n = np.zeros(nbins)
	nmt = np.zeros(nbins)
	if drag:
		nd = np.zeros(nbins)
	if triaxial:
		nt = np.zeros(nbins)
	
	nsubs = 0
	ri = 0
	for i in range(len(subs)-1):
		# if host.subhalos[i]:
		# print host.subhalos[i]
		zz = 1/subs[i].a - 1
		tt = host.cosmo.age(zz)
		vmax = subs[i].v_max[-1]
		if vmax > v_thresh and len(tt) > 10:
			hostM = host.M_sp(tt)
			hostR = host.virial_radius(hostM, zz)
			dd = np.sqrt(subs[i].rel_x*subs[i].rel_x + subs[i].rel_y*subs[i].rel_y + subs[i].rel_z*subs[i].rel_z) # Mpc
			ratio_to_tt = interp1d(dd/hostR,tt)

			try:
				t0 = ratio_to_tt(2) # determine when subhalo is at a distance of 2 * host.R
			except:
				continue

			nsubs += 1

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
			if triaxial:
				diff = np.abs(rt[ri] - rbc*host.R)
				rbin = diff.argmin()
				nt[rbin] += 1	
			ri+=1
	print 'nsubs = %i' %nsubs

	vshell = 4/3*np.pi*((rbins[1:]*host.R)**3-(rbins[:-1]*host.R)**3)
	
	sd = n*dr/vshell
	sd_mt = nmt*dr/vshell
	if drag: sd_d = nd*dr/vshell
	if triaxial: sd_t = nt*dr/vshell

	sigma[j] = sd
	sigma_mt[j] = sd_mt
	if drag: sigma_d[j] = sd_d
	if triaxial: sigma_t[j] = sd_t

	# sigma*=sd 
	# sigma_mt*=sd_mt
	# sigma_d*=sd_d

	# var+=(nsubs/sd**2)
	# var_mt+=(nsubs/sd_mt**2)
	# var_d+=(nsubs/sd_d**2)

sigma = np.median(sigma,axis=0)
sigma_mt = np.median(sigma_mt,axis=0)
if drag: sigma_d = np.median(sigma_d,axis=0)
if triaxial: sigma_t = np.median(sigma_t,axis=0)

# print 'results'
# print sigma

# err = sigma/nhosts*np.sqrt(var)
# err_mt = sigma/nhosts*np.sqrt(var_mt)
# err_d = sigma/nhosts*np.sqrt(var_d)

# PLOTTING
fig, ax = plt.subplots(1,2,figsize=(20,8))
error = 0
if error:
	ax[0].errorbar(rbc,sigma,xerr=dr/2,yerr=err,ls='-',color='b',label=r'$\mathrm{Spherical\ NFW}$',lw=3,markersize='10')
	ax[0].errorbar(rbc,sigma_d,xerr=dr/2,yerr=err_d,ls='-',color='r',label=r'$\mathrm{Drag\ Force}$',lw=3,markersize='10')
	ax[0].errorbar(rbc,sigma_mt,xerr=dr/2,yerr=err_mt,ls='-',color='g',label=r'$\mathrm{Merger\ Tree}$',lw=3,markersize='10')
else:
	ax[0].plot(rbc,sigma,'-',color='r',label=r'$\mathrm{Spherical\ NFW}$',lw=3,markersize='10')
	if drag: ax[0].plot(rbc,sigma_d,'-',color='b',label=r'$\mathrm{Drag\ Force}$',lw=3,markersize='10')
	if triaxial: ax[0].plot(rbc,sigma_t,'-',color='m',label=r'$\mathrm{Triaxial\ NFW}$',lw=3,markersize='10')
	ax[0].plot(rbc,sigma_mt,'-',color='g',label=r'$\mathrm{Merger\ Tree}$',lw=3,markersize='10')
	
ax[1].plot(rbc,sigma_d/sigma,'-',color='b',label=r'$\mathrm{Drag\ Force\ Ratio}$',lw=3)

ax[0].grid()
ax[0].set_xlabel(r'$\mathrm{r/R_{200m}}$')
ax[0].set_ylabel(r'$\mathrm{\Sigma /R_{200m}\ (subhalos/Mpc^3)}$')
ax[0].set_title(r'$\mathrm{%i\ Hosts\ Stacked,\ v_{thresh}\ =\ %.2f\ km/s}$' %(nhosts,v_thresh))
ax[0].legend()
ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].set_xlim(0,3*host.R)

ax[1].grid()
ax[1].set_xlabel(r'$\mathrm{r/R_{200m}}$')
ax[1].set_ylabel(r'$\mathrm{\Sigma_{drag} /\Sigma}$')
ax[1].legend()
# ax[1].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_xlim(0,3*host.R)

plt.savefig('plots/subhalo_distribution_%.0e_%i_%i.png'  %(dt,v_thresh,nhosts))

# plt.show()
