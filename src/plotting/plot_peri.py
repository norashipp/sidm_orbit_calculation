#!/usr/bin/env python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy import loadtxt
from scipy.signal import savgol_filter

dpi = 175
fontsize = 15
plt.rc('savefig', dpi=dpi)
plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.rc('xtick.major', pad=5)
plt.rc('xtick.minor', pad=5)
plt.rc('ytick.major', pad=5)
plt.rc('ytick.minor', pad=5)

opt = 'apocenter'
opti = 'a'
# opt = 'pericenter'
# opti = 'p'

dt = 4e-3
potential = 'spherical_NFW'
nhosts = 51

# bmin = -0.05
bmax = 1.1
# bins = np.linspace(0.0,0.8,20)
bins = np.logspace(np.log10(0.0001),np.log10(bmax),20)
# bins = 10

colors = np.array(['b','g','r','c','y','m','k','orange'],dtype=str)

for i,sig in enumerate([3,9,15,21]):
	c = colors[i]
	dperi = loadtxt('output/first_%s_%s_%.0e_sigma_%.2f.txt' %(opt,potential,dt,sig))
	dperi = dperi[dperi>0]
	large = dperi[dperi>=0.1]
	small = dperi[dperi<0.1]
	print 'sigma = %i' %sig
	print 'large = ', len(large)
	print 'small = ', len(small)
	print 'fraction = ', len(large)/len(small)
	print 'median = ', np.median(dperi)
	print 
	# dperi = savgol_filter(dperi,15,2)
	# dperi = savgol_filter(np.log10(dperi),15,2)
	# dperi = 10**dperi
	plt.hist(dperi,color=c,bins=bins,histtype='step',lw=3,label=r'$\sigma/m_{\rm \chi} = %i\ \mathrm{cm^2/g}$'%sig)
	plt.plot([np.median(dperi),np.median(dperi)],[0,1.4],'-',color=c,lw=3)
	# plt.plot(np.median(dperi),0,'+',color=c,lw=3,ms=20)

plt.xlabel(r'$\Delta r_{\rm %s} / R_{\rm 200m}$' %opti)
plt.ylabel(r'$N_{\rm subhalos}$')
plt.title(r'$\mathrm{Change\ in\ First\ Apocenter\ of\ Subhalo\ Orbits}$')
plt.yscale('log')
plt.xscale('log')
# plt.ylim(0,1e2)
# plt.legend()
plt.grid()
plt.savefig('plots/first_%s_%s_%.0e_%i.png'  %(opt,potential,dt,nhosts))
plt.show()
