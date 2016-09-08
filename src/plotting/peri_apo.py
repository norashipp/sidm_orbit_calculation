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
	if len(minima) < 2 or len(maxima) < 2:
		return None
	else:
		if r[1] > r[0]:
			rp = minima[1]
			ra = maxima[0]
		elif r[1] < r[0]:
			rp = minima[0]
			ra = maxima[1]
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

rperi = [] 
rapo = []
rperi_mt = [] 
rapo_mt = []
rperi_d = [] 
rapo_d = []

ratio = []
ratio_mt = []
ratio_d = []

for j in range(nhosts):
	host_idx = hosts[j]
	host = HostHalo(host_idx,potential)
	host.update(host.cosmo.age(0))
	print 'Host %i' %host_idx
	print 'M = %.2e' %host.M
	subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)
	
	for i in range(len(host.subhalos)):
		print 'subhalo %i' %i
		if host.subhalos[i]:
			sub_idx = i
			print 'calculating subhalo %i' %i
			# infile = HOMEDIR+'sidm_orbit_calculation/src/output/leapfrog_%s/%i_leapfrog_%s_%.0e_%i.dat' %(potential,host_idx,potential,dt,sub_idx)
			infile = HOMEDIR+'sidm_orbit_calculation/src/output/%i_leapfrog_%s_%.0e_%i.dat' %(host_idx,potential,dt,sub_idx)
			f = open(infile,'rb')
			data = cPickle.load(f)
			f.close()
			t, positions, _ = data
			d = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)

			if drag:
				sigma = 6
				# infile = HOMEDIR+'sidm_orbit_calculation/src/output/dissipative_%s/sigma_%i/%i_dissipative_%s_%.0e_%i.dat' %(potential,sigma,host_idx,potential,dt,sub_idx)
				infile = HOMEDIR+'sidm_orbit_calculation/src/output/sigma%i/%i_dissipative_%s_%.0e_%i.dat' %(sigma,host_idx,potential,dt,sub_idx)
				f = open(infile,'rb')
				data = cPickle.load(f)
				f.close()
				td, positions, _ = data
				dd = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)

			sub = host.subhalos[i]
			zz = 1/subs[i].a - 1
			tt = host.cosmo.age(zz)
			vmax = subs[i].v_max[-1]
			if vmax > v_thresh:
				# merger tree
				dmt = np.sqrt(subs[i].rel_x**2 + subs[i].rel_y**2 + subs[i].rel_z**2)
				dmt = dmt[tt>sub.t0]
				tt = tt[tt>sub.t0]
				
				rpmt,ramt = get_radii(dmt)
				if rpmt == None: continue
				ratio_mt.append(ramt/rpmt)

				# spherical NFW
				rp,ra = get_radii(d)
				ratio.append(ra/rp)

				# drag force
				if drag:
					rpd,rad = get_radii(dd)
					ratio.append(rad/rpd)

				plt.figure()
				plt.plot(tt,dmt)
				plt.plot([tt.min(),tt.max()],[rpmt,rpmt],'*',markersize=10)
				plt.plot([tt.min(),tt.max()],[ramt,ramt],'*',markersize=10)
				plt.show()
				break
'''
rperi = np.asarray(rperi)
rapo = np.asarray(rapo)
rperi_mt = np.asarray(rperi_mt)
rapo_mt = np.asarray(rapo_mt)
rperi_d = np.asarray(rperi_d)
rapo_d = np.asarray(rapo_d)
'''

ratio = np.asarray(ratio)
ratio_mt = np.asarray(ratio_mt)
ratio_d = np.asarray(ratio_d)

print ratio.min(), ratio.max()
print ratio_d.min(), ratio_d.max()

# PLOTTING
# bins = np.linspace(0,2,50)
# bins=np.linspace(1,3,20)
plt.figure()
plt.hist(ratio,histtype='step',lw=3,normed=True,label='spherical NFW')
# plt.hist(ratio_mt,bins=bins,lw=3,normed=True,label='merger tree',alpha=0.3)
plt.hist(ratio_d,histtype='step',lw=3,normed=True,label='drag')
plt.grid()
plt.xlabel(r'$\mathrm{r/R_{200m}}$')
plt.ylabel(r'$\mathrm{subhalos}$')
# plt.title(r'$\mathrm{%i\ Hosts\ Stacked,\ v_{thresh}\ =\ %.2f\ km/s}$' %(nhosts,v_thresh))
plt.legend()
# plt.yscale('log')
plt.savefig('plots/pericenter_apocenter_ratio_%s_%.0e_%i.png'  %(potential,dt,nhosts))
plt.show()
