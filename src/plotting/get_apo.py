#!/usr/bin/env python

from __future__ import division
import numpy as np
import sys
import glob
import cPickle

from sidm_orbit_calculation.src.utils.setup import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.merger_tree.cluster import *

# first after infall
def get_radii(r):
	minima = np.r_[True, r[1:] < r[:-1]] & np.r_[r[:-1] < r[1:], True]
	maxima = np.r_[True, r[1:] > r[:-1]] & np.r_[r[:-1] > r[1:], True]
	if len(r[minima]) < 2 or len(r[maxima]) < 2:
		return -1, -1
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
		return -1, -1
	else:
		if r[-1] < r[-2]:
			rp = r[minima][-2]
			ra = r[maxima][-1]
		elif r[-1] > r[-2]:
			rp = r[minima][-1]
			ra = r[maxima][-2]
		return rp, ra

hosts = np.array(sys.argv[1:],dtype=int)
if len(hosts) == 0: hosts = np.arange(51)
nhosts = len(hosts)
print '%i hosts' %nhosts

potential = 'spherical_NFW'
dt = 4e-3

host_R = np.loadtxt(HOMEDIR+'sidm_orbit_calculation/src/plotting/host_R.txt')

sigs = [0,3,9,15,21]
for sig in sigs:
	apo = []
	peri = []
	for host_idx in hosts:
                print 'Host %i' %host_idx
		R = host_R[host_idx]
		infiles = glob.glob(HOMEDIR+'data/candidacy/sigma%i/%i_%s_%.0e_%.2f_*.dat' %(sig,host_idx,potential,dt,sig))	
		for infile in infiles:
			f = open(infile,'rb')
			data = cPickle.load(f)
			f.close()
			_, positions, _ = data
			d = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)/R
			rp,ra = get_radii(d)
			apo.append(ra)
			peri.append(rp)
			
	np.savetxt('output/apocenter_%s_%.0e_sigma_%.2f_nonorm.txt' %(potential, dt, sig),apo)
	np.savetxt('output/pericenter_%s_%.0e_sigma_%.2f_nonorm.txt' %(potential, dt, sig),peri)
