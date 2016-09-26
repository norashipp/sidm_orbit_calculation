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
	if len(r[minima]) < 2 or len(r[maxima]) < 2:
		return None, None
	else:
		if r[1] < r[0]:
			rp = r[minima][1]
			ra = r[maxima][0]
		elif r[1] > r[0]:
			rp = r[minima][0]
			ra = r[maxima][1]
		return rp, ra
# last
def get_final_radii(r):
	minima = np.r_[True, r[1:] < r[:-1]] & np.r_[r[:-1] < r[1:], True]
	maxima = np.r_[True, r[1:] > r[:-1]] & np.r_[r[:-1] > r[1:], True]
	if len(r[minima]) < 2 or len(r[maxima]) < 2:
		return None, None
	else:
		if r[-1] < r[-2]:
			rp = r[minima][-2]
			ra = r[maxima][-1]
		elif r[-1] > r[-2]:
			rp = r[minima][-1]
			ra = r[maxima][-2]
		return rp, ra

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

apo0 = []
peri0 = []
for j, host_idx in enumerate(hosts):
	host = HostHalo(host_idx,potential,subs=True,scale_density=False)
        host.update(host.cosmo.age(0))
	subhalos = np.copy(host.subhalos)
	sub_dict[j] = subhalos
	print 'Host %i' %host_idx
	print 'M = %.2e' %host.M
        
        for sub_idx, sub in enumerate(subhalos):
            if sub:
                infile = HOMEDIR+'data/candidacy/sigma%i/%i_%s_%.0e_%.2f_%i.dat' %(0,host_idx,potential,dt,0.,sub_idx)
                f = open(infile,'rb')
                data = cPickle.load(f)
                f.close()
                _, positions, _ = data
                d = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)/host.R
                rp0,ra0 = get_radii(d)
		if rp0 == None or ra0 == None: continue
                apo0.append(ra0)
                peri0.append(rp0)

a95 = np.percentile(apo0,95)

sigs = [3,9,15,21]
for sig in sigs:
	dapo = []
	dperi = []
	for j, host_idx in enumerate(hosts):
		host = HostHalo(host_idx,potential,subs=False,scale_density=False)
		subhalos = sub_dict[j]
		host.update(host.cosmo.age(0))
		for sub_idx, sub in enumerate(subhalos):
			if sub:
				infile = HOMEDIR+'data/candidacy/sigma%i/%i_%s_%.0e_%.2f_%i.dat' %(0,host_idx,potential,dt,0.,sub_idx)
				f = open(infile,'rb')
				data = cPickle.load(f)
				f.close()
				_, positions, _ = data
				d = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)
				rp0,ra0 = get_radii(d)
				
				infile = HOMEDIR+'data/candidacy/sigma%i/%i_%s_%.0e_%.2f_%i.dat' %(sig,host_idx,potential,dt,sig,sub_idx)
				f = open(infile,'rb')
				data = cPickle.load(f)
				f.close()
				_, positions, _ = data
				d = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)/host.R
				rp,ra = get_radii(d)
				if ra0 == None or ra ==None or ra0 < a95: continue
				dapo.append(ra0-ra)
				dperi.append(rp0-rp)

	np.savetxt('output/apocenter_a95_%s_%.0e_sigma_%.2f.txt' %(potential, dt, sig),dapo)
	np.savetxt('output/pericenter_a95_%s_%.0e_sigma_%.2f.txt' %(potential, dt, sig),dperi)
