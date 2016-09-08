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

hosts = np.array(sys.argv[1:],dtype=int)
if len(hosts) == 0: hosts = np.arange(51)
nhosts = len(hosts)
print '%i hosts' %nhosts

potential = 'spherical_NFW'
dt = 4e-3
drag = 1
v_thresh = 100 # km/s

nbins = 12
rbins = np.linspace(0,3,nbins+1)
# dr = rbins[1]-rbins[0]
# rbc = rbins[1:]-dr/2


sigs = [3,6]
colors = ['b','g','r','c','m','y','k','w']
for j, sig in enumerate(sigs):
	dperi = []
	dapo = []
	
	for j in range(nhosts):
		host_idx = hosts[j]
		host = HostHalo(host_idx,potential)
		host.update(host.cosmo.age(0))
		print 'Host %i' %host_idx
		print 'M = %.2e' %host.M
		subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)
		
		plt.figure()

		for i in range(len(host.subhalos)):
			# print 'subhalo %i' %i
			if host.subhalos[i]:
				sub_idx = i
				# print 'calculating subhalo %i' %i
				# infile = HOMEDIR+'sidm_orbit_calculation/src/output/leapfrog_%s/%i_leapfrog_%s_%.0e_%i.dat' %(potential,host_idx,potential,dt,sub_idx)
				infile = HOMEDIR+'sidm_orbit_calculation/src/output/leapfrog_%s/%i_leapfrog_%s_%.0e_%i_major_axis.dat' %(potential,host_idx,potential,dt,sub_idx)
				f = open(infile,'rb')
				data = cPickle.load(f)
				f.close()
				t, positions, _ = data
				d = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)/host.R
				# spherical NFW
				rp,ra = get_radii(d)
				if not rp: continue
				# rapo.append(ra)
				# rperi.append(rp)

				'''
				# merger tree
				dmt = np.sqrt(subs[i].rel_x**2 + subs[i].rel_y**2 + subs[i].rel_z**2)/host.R
				dmt = dmt[tt>sub.t0]
				tt = tt[tt>sub.t0]
				rpmt,ramt = get_radii(dmt)
				if rpmt == None: continue
				rapo_mt.append(ramt)
				rper_mt.append(rpmt)

				plt.figure()
				plt.plot(tt,dmt)
				plt.plot([tt.min(),tt.max()],[rpmt,rpmt],'*',markersize=10)
				plt.plot([tt.min(),tt.max()],[ramt,ramt],'*',markersize=10)
				plt.show()
				break
				'''
			
				# drag force
				# infile = HOMEDIR+'sidm_orbit_calculation/src/output/dissipative_%s/sigma_%i/%i_dissipative_%s_%.0e_%i.dat' %(potential,sigma,host_idx,potential,dt,sub_idx)
				infile = HOMEDIR+'sidm_orbit_calculation/src/output/dissipative_%s/sigma_%i/%i_dissipative_%s_%.0e_%i_major_axis.dat' %(potential,sig,host_idx,potential,dt,sub_idx)
				f = open(infile,'rb')
				data = cPickle.load(f)
				f.close()
				td, positions, _ = data
				dd = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)/host.R

				rpd,rad = get_radii(dd)
				# rapo_d.append(rad)
				# rperi_d.append(rperi)

				dp = rp-rpd
				da = ra-rad

				# print rp, rpd, dp
				# print ra, rad, da
				# print 

				dperi.append(dp)
				dapo.append(da)

		c = colors[j]
		# FIX COLORS
		# print dperi
		# print dapo
		plt.hist(dperi,histtype='step',lw=3,normed=True,label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sig)
		# plt.hist(dapo,histtype='step',lw=3,normed=True)

	plt.grid()
	plt.xlabel(r'$\mathrm{\Delta r_{p}/R_{200m}}$')
	# plt.xlabel(r'$\mathrm{\Delta r_{p}/R_{200m}}$')
	plt.ylabel(r'$\mathrm{subhalos}$')
	plt.legend()
	# plt.savefig('plots/pericenter_%s_%.0e_%i.png'  %(potential,dt,nhosts))

plt.show()
