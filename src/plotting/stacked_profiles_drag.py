from __future__ import division
import pylab
import numpy as np
import cPickle
import sys
import glob
from scipy.interpolate import interp1d
import time
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.utils.setup import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.merger_tree.cluster import *

start_time = time.time()

dpi = 175
fontsize = 20
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

potential = 'spherical_NFW'
dt = 4e-3

mt = 0

v_thresh = 70 # km/s

nbins = 18
rbins = np.logspace(np.log10(0.1),np.log10(3),nbins+1)
dr = rbins[1:]-rbins[:-1]
rbc = rbins[1:]-dr/2

rbc_plt = np.logspace(np.log10(0.1),np.log10(3),nbins*3)

sub_dict = {}

sigs = [0,3,9,15,21]
# sigs = [0,3,6,9,12,15,18,21]
colors = np.array(['b','g','r','c','y','m','k','orange'],dtype=str)

fig, ax = plt.subplots(1,2,figsize=(20,8))
	
for i, sig in enumerate(sigs):
	sigma = np.ones((nhosts,nbins))
	if mt: sigma_mt = np.ones((nhosts,nbins))
	c = colors[i]
	for j, host_idx in enumerate(hosts):
		if sig == 0:
			host = HostHalo(host_idx,potential,subs=True,scale_density=False)
			subhalos = np.copy(host.subhalos)
			sub_dict[j] = subhalos
		else:
			host = HostHalo(host_idx,potential,subs=False,scale_density=False)
			subhalos = sub_dict[j]

		host.update(host.cosmo.age(0))
		
		if mt:
			subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)
			nmt = np.zeros(nbins)
			for i in range(len(subs)-1):
				rf = np.sqrt(subs[i].rel_x[-1]**2 + subs[i].rel_y[-1]**2 + subs[i].rel_z[-1]**2)
				diff = np.abs(rf - rbc*host.R)
				rbin = diff.argmin()
				nmt[rbin] += 1

		infile = HOMEDIR+'Dropbox/SIDMdata/candidacy/final_positions/sigma%i/final_positions_%i_%s_%.0e_%.2f.txt' %(sig,host_idx,potential,dt,sig)
		x,y,z = np.loadtxt(infile,unpack=True)
		r = np.sqrt(x**2 + y**2 + z**2)

		print 'Sigma = %i' %sig
		print 'Host %i' %host_idx
		print 'M = %.2e' %host.M
		print 'Number of subhalos = %i' %len(r)

		n = np.zeros(nbins)		
		nsubs = 0
		ri = 0
		
		for sub in subhalos:
			if sub:
				diff = np.abs(r[ri] - rbc*host.R)
				rbin = diff.argmin()
				n[rbin] += 1
				nsubs += 1
				ri+=1

		print 'nsubs = %i' %nsubs
		print 

		vshell = 4/3*np.pi*((rbins[1:]*host.R)**3-(rbins[:-1]*host.R)**3)
		
		sd = n*dr/vshell
		if mt: sd_mt = n_mt*dr/vshell

		sigma[j] = sd
		if mt: sigma_mt[j] = sd_mt

	sigma = np.median(sigma,axis=0)

	F = open('stacked_profile_nbins_18_%.2f.txt' %sig, 'wb')
	cPickle.dump(sigma,F)
	# np.savetxt(fname,sigma)

	sigma_plt = savgol_filter(sigma,3,2)

	if mt: sigma_mt = np.median(sigma_mt,axis=0) # add spline if using this
	if sig == 0: sigma0_plt = np.copy(sigma_plt)

'''
	# PLOTTING
	ax[0].plot(rbc,sigma_plt,'-',lw=3,color=c,label=r'$\sigma/m_{\chi} = %i\ \mathrm{cm^2/g}$' %sig)		
	if sig > 0:
		ax[1].plot(rbc,sigma_plt/sigma0_plt,'-',lw=3,color=c,label=r'$\sigma/m_{\chi} = %i\ \mathrm{cm^2/g}$' %sig)


if mt: ax[0].plot(rbc,sigma_mt,'-',color='g',label=r'$\mathrm{Merger\ Tree}$',lw=3,markersize='10')

ax[0].set_xlabel(r'$r/R_{\rm 200m}$')
ax[0].set_ylabel(r'$\Sigma /R_{\rm 200m}\ \mathrm{(subhalos/Mpc^3)}$')
ax[0].set_title(r'$\mathrm{%i\ Hosts\ Stacked,\ v_{thresh}\ =\ %.2f\ km/s}$' %(nhosts,v_thresh))
# ax[0].legend(fontsize=18)
ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].set_xlim(0,3)
ax[0].grid()

ax[1].set_xlabel(r'$r/R_{\rm 200m}$')
ax[1].set_ylabel(r'$\Sigma_{\rm drag} / \Sigma$')
# ax[1].legend(fontsize=18)
# ax[1].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_xlim(0,3)
ax[1].grid()

plt.savefig('plots/subhalo_distribution_%.0e_%i_%i.png'  %(dt,v_thresh,nhosts))

end_time = time.time()
print 'Time passed = %.2f' %(end_time-start_time)

# plt.show()
'''