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

'''
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
# rbins = np.linspace(0,3,nbins+1)
rbins = np.logspace(-4,-1,nbins+1)
# dr = rbins[1]-rbins[0]
# rbc = rbins[1:]-dr/2

plt.figure()

sigs = [3,9,15,21]
colors = ['b','g','r','c','m','y','k','w']

'''
for j in range(nhosts):
	host_idx = hosts[j]
	host = HostHalo(host_idx,potential,subs=True,scale_density=False)
        host.update(host.cosmo.age(0))
	subhalos = np.copy(host.subhalos)
	sub_dict[j] = subhalos

	print 'Host %i' %host_idx
	print 'M = %.2e' %host.M

	infile = HOMEDIR+'data/candidacy/sigma0/%i_%s_%.0e_0.00_%i.dat' %(host_idx,potential,dt,sub_idx)
	f = open(infile,'rb')
	data = cPickle.load(f)
	f.close()
	_, positions, _ = data
	d = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)/host.R
	rp,ra = get_radii(d)
	if not rp: continue
'''

for k, sig in enumerate(sigs):
	dperi = []
	dapo = []
	
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
                infile = HOMEDIR+'data/candidacy/sigma0/%i_%s_%.0e_0.00_%i.dat' %(host_idx,potential,dt,sub_idx)
                f = open(infile,'rb')
                data = cPickle.load(f)
                f.close()
                _, positions, _ = data
                d = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)/host.R
                rp,ra = get_radii(d)
                if not rp: continue
				
                infile = HOMEDIR+'data/candidacy/sigma%i/%i_%s_%.0e_%.2f_%i.dat' %(sig,host_idx,potential,dt,sig,sub_idx)
				f = open(infile,'rb')
				data = cPickle.load(f)
				f.close()
				td, positions, _ = data
				dd = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)/host.R

				rpd,rad = get_radii(dd)
                if not rpd: continue
                if rpd > 0.2: continue
				dp = rp-rpd
				da = ra-rad

				dperi.append(dp)
				dapo.append(da)
	
	writing = 1
	if writing:
		np.savetxt('output/pericenter_%s_%.0e_sigma_%.2f.txt' %(potential, dt, sig),dperi)
		np.savetxt('output/apocenter_%s_%.0e_sigma_%.2f.txt' %(potential, dt, sig),dapo)
	
    print rbins
	print min(dperi), max(dperi)
        plt.hist(dperi,bins=rbins,histtype='step',lw=3,label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sig)

	plt.xlabel(r'$\mathrm{\Delta r_{p}/R_{200m}}$')
	plt.ylabel(r'$\mathrm{subhalos}$')
	plt.yscale('log')
    plt.xscale('log')
	plt.legend()
	plt.grid()
	plt.savefig('plots/pericenter_%s_%.0e_%i.png'  %(potential,dt,nhosts))

# plt.show()
