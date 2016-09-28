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

dt = 4e-3
potential = 'spherical_NFW'
nhosts = 51

bmax = 0.5
bins = np.logspace(np.log10(0.001),np.log10(bmax),20)
# bins = np.logspace(np.log10(0.0001),np.log10(bmax),20)

colors = np.array(['b','g','r','c','y','m','k','orange'],dtype=str)

apo0 = loadtxt('output/%s_%s_%.0e_sigma_0.00_nonorm.txt' %(opt,potential,dt))
p95 = np.percentile(apo0[apo0 > 0],95)
p75 = np.percentile(apo0[apo0 > 0],75)
# print p95
# idx = (apo0 >= p95)
# apo95 = apo0[idx]
# print len(apo95)
# print len(idx)
# print len(apo0)

for i,sig in enumerate([3,9,15,21]):
	c = colors[i]
	apo = loadtxt('output/%s_%s_%.0e_sigma_%.2f_nonorm.txt' %(opt,potential,dt,sig))
	# print len(apo), len(idx)
	# idx = ((apo0 >= p75) & (apo > 0))
	idx = ((apo0 > 0) & (apo > 0))
	
	apo = apo[idx]
	apo95 = apo0[idx]
	

	print 'sigma = %i - %i subhalos' %(sig,len(apo))
	
	dapo = apo95-apo
	dapo = dapo[dapo>0]
	
	dapo.sort()
	# print dapo

	print dapo.min(), dapo.max()
	
	meds = []
	for i in range(5000):
		a = np.random.choice(dapo,size=len(dapo))
		meds.append(np.median(a))
	std = np.std(meds)
	print 'standard deviation = ', std
	print 'mean = ', np.mean(meds)
	print 

	plt.hist(dapo,color=c,bins=bins,histtype='step',lw=3,label=r'$\sigma/m_{\rm \chi} = %i\ \mathrm{cm^2/g}$'%sig)
	plt.plot([np.median(dapo),np.median(dapo)],[0,1.4],'-',color=c,lw=3)

plt.xlabel(r'$\Delta r_{\rm %s} / R_{\rm 200 m}$' %(opti))
plt.ylabel(r'$N_{\rm subhalos}$')
plt.title(r'$\mathrm{Change\ in\ First\ Apocenter\ of\ Subhalo\ Orbits}$')
plt.yscale('log')
plt.xscale('log')
# plt.ylim(0,1e2)
plt.legend()
plt.grid()
plt.savefig('plots/first_%s_fraction_%s_%.0e_%i.png'  %(opt,potential,dt,nhosts))
plt.show()
