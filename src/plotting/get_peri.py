#!/usr/bin/env python

from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import cPickle

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

# first after infall
def get_radii(r):
	minima = np.r_[True, r[1:] < r[:-1]] & np.r_[r[:-1] < r[1:], True]
	maxima = np.r_[True, r[1:] > r[:-1]] & np.r_[r[:-1] > r[1:], True]
	# print r[minima]
	# print r[maxima]
	if len(r[minima]) < 2 or len(r[maxima]) < 2:
		# print 'problem..'
		return None, None
	else:
		if r[1] < r[0]:
			rp = r[minima][1]
			ra = r[maxima][0]
		elif r[1] > r[0]:
			rp = r[minima][0]
			ra = r[maxima][1]
		return rp, ra
'''
# last
def get_radii(r):
	minima = np.r_[True, r[1:] < r[:-1]] & np.r_[r[:-1] < r[1:], True]
	maxima = np.r_[True, r[1:] > r[:-1]] & np.r_[r[:-1] > r[1:], True]
	# print r[minima]
	# print r[maxima]
	if len(r[minima]) < 2 or len(r[maxima]) < 2:
		# print 'problem..'
		return None, None
	else:
		if r[-1] < r[-2]:
			rp = r[minima][-2]
			ra = r[maxima][-1]
		elif r[-1] > r[-2]:
			rp = r[minima][-1]
			ra = r[maxima][-2]
		return rp, ra
'''

sub_dict = {}
hosts = np.array(sys.argv[1:],dtype=int)
if len(hosts) == 0: hosts = np.arange(51)
nhosts = len(hosts)
print '%i hosts' %nhosts

potential = 'spherical_NFW'
dt = 4e-3
drag = 1
v_thresh = 100 # km/s

nbins = 12
rbins = np.logspace(-4,-1,nbins+1)

sigs = [0,3,9,15,21]
colors = ['b','g','r','c','m','y','k','w']
for k, sig in enumerate(sigs):
	apo = []
	peri = []
	for j, host_idx in enumerate(hosts):
		if sig == sigs[0]:
			host = HostHalo(host_idx,potential,subs=True,scale_density=False)
			subhalos = np.copy(host.subhalos)
			sub_dict[j] = subhalos
		else:
			host = HostHalo(host_idx,potential,subs=False,scale_density=False)
			subhalos = sub_dict[j]
		host.update(host.cosmo.age(0))
		print 'Host %i' %host_idx
		print 'M = %.2e' %host.M
		
		for sub_idx, sub in enumerate(subhalos):
			if sub:
				infile = HOMEDIR+'data/candidacy/sigma%i/%i_%s_%.0e_%.2f_%i.dat' %(sig,host_idx,potential,dt,sig,sub_idx)
				f = open(infile,'rb')
				data = cPickle.load(f)
				f.close()
				_, positions, _ = data
				d = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)
				rp,ra = get_radii(d)
				# if not rp: continue
				apo.append(ra)
				peri.append(rp)

	np.savetxt('output/first_apocenter_%s_%.0e_sigma_%.2f.txt' %(potential, dt, sig),apo)
	np.savetxt('output/first_pericenter_%s_%.0e_sigma_%.2f.txt' %(potential, dt, sig),peri)
