#!/usr/bin/env python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy import loadtxt

dpi = 175
fontsize = 15
plt.rc('savefig', dpi=dpi)
plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.rc('xtick.major', pad=5)
plt.rc('xtick.minor', pad=5)
plt.rc('ytick.major', pad=5)
plt.rc('ytick.minor', pad=5)

# opt = 'apocenter'
# opti = 'a'
opt = 'pericenter'
opti = 'p'

dt = 4e-3
potential = 'spherical_NFW'
nhosts = 51

bmin = -0.05
bmax = 0.12
bins = np.linspace(-0.05,0.12,20)
# bins = np.logspace(np.log10(0.0001),np.log10(bmax),20)
bins = 10
colors = np.array(['b','g','r','c','y','m','k','orange'],dtype=str)
for i,sig in enumerate([3,9,15,21]):
	c = colors[i]
	dperi = loadtxt('output/first_%s_%s_%.0e_sigma_%.2f.txt' %(opt,potential,dt,sig))
	plt.hist(dperi,color=c,bins=bins,histtype='step',lw=3,label=r'$\sigma/m_{\rm \chi} = %i\ \mathrm{cm^2/g}$'%sig)
	plt.plot([np.median(dperi),np.median(dperi)],[0,1e3],'--',color=c,lw=3)

plt.xlabel(r'$\Delta r_{\rm %s} / R_{\rm 200m}$' %opti)
plt.ylabel(r'$N_{\rm subhalos}$')
plt.title(r'$\mathrm{Change\ in\ First\ Pericenter\ of\ Subhalo\ Orbits}$')
plt.yscale('log')
# plt.xscale('log')
plt.legend()
plt.grid()
plt.savefig('plots/first_%s_%s_%.0e_%i.png'  %(opt,potential,dt,nhosts))
plt.show()
