import cPickle
import numpy as np
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import sys
from numpy import loadtxt

from sidm_orbit_calculation.src.utils.setup import *
from sidm_orbit_calculation.src.halos.host_halo import *
import sidm_orbit_calculation.src.merger_tree.cluster as mt

dpi = 175
fontsize = 15
plt.rc('savefig', dpi=dpi)
plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.rc('xtick.major', pad=5)
plt.rc('xtick.minor', pad=5)
plt.rc('ytick.major', pad=5)
plt.rc('ytick.minor', pad=5)

def kde_sklearn(dista, grid, bandwidth = 1.0, **kwargs):
    kde_skl = KernelDensity(bandwidth = bandwidth, **kwargs)
    kde_skl.fit(dista[:, np.newaxis])
    log_pdf = kde_skl.score_samples(grid[:, np.newaxis]) # sklearn returns log(density)

    return np.exp(log_pdf)

hosts = np.array(sys.argv[1:],dtype=int)
integrator = 'leapfrog'
dt = 4e-3

for host_idx in hosts:
	potential = 'spherical_NFW'
	# dist = loadtxt(HOMEDIR + 'sidm_orbit_calculation/src/output/final_distances_%i.txt' %host_idx, unpack=True)
	pos = loadtxt(HOMEDIR+'sidm_orbit_calculation/src/output/final_positions_%i_%s_%s_%.0e.txt' %(host_idx,integrator,potential,dt))
	dist = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
	
	host = HostHalo(host_idx,potential)
	host.update(host.cosmo.age(0))
	host_radius = host.R

	triaxial = 0
	if triaxial:
		potential = 'triaxial_NFW'
		pos = loadtxt(HOMEDIR+'sidm_orbit_calculation/src/output/final_positions_%i_%s_%s_%.0e.txt' %(host_idx,integrator,potential,dt))
		dist_tri = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
		print dist.shape, dist_tri.shape, pos.shape

	compare = 1
	if compare:
		subs = mt.SubHalos(HOMEDIR + 'sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat' % host_idx)
		dist_mt = []
		mass = []
		for i in range(len(subs)-1):
			sub = subs[i]
			r = np.sqrt(sub.rel_x[-1]**2 + sub.rel_y[-1]**2 + sub.rel_z[-1]**2)
			dist_mt.append(r)
			mass.append(sub.m_200m[-1])
		dist_mt = np.asarray(dist_mt)[np.not_equal(host.subhalos,None)]
		mass = np.asarray(mass)[np.not_equal(host.subhalos,None)]
		# print dist_mt.shape, mass.shape, dist.shape

	plot_kde = 0
	if plot_kde:
		sigma = np.std(dist)
		sigma_mt = np.std(dist_mt)
		print 'std = ', sigma, sigma_mt
		n = len(dist)
		h = 1.059*sigma*n**-(1/5)
		h_mt = 1.059*sigma_mt*n**-(1/5)
		print 'h = ', h
		print 'h_mt = ', h_mt

		grid = np.arange(0,20,0.1)
		kde = kde_sklearn(dist/host_radius, grid, h, kernel='gaussian')
		kde_mt = kde_sklearn(dist_mt/host_radius, grid, h_mt, kernel='gaussian')

		plt.figure()
		plt.plot(grid,kde,'c',lw=3)
		plt.plot(grid,kde_mt,'g',lw=3)
		plt.title(r'$\mathrm{Host\ %i\ -\ Gaussian\ KDE}$' %host_idx)
		plt.xlabel(r'$\mathrm{r/R_{200}}$')
		plt.ylabel(r'$\mathrm{subhalos}$')
		plt.grid()
		plt.xlim([0,10])
		# plt.savefig('%i_distribution_kde.png' %host_idx)
	
	# q75, q25 = np.percentile(dist, [75 ,25])
	# iqr = q75 - q25
	# h = 2*iqr*len(dist)**(-1/3)
	# nbins = (dist.max()-dist.min())/h
	# nbins = np.sqrt(len(dist))
	# print nbins
	nbins = 10

	plt.figure()
	plt.hist(dist,color='c',bins=nbins,histtype='step',normed=True,lw=3,weights=mass,label=r'$\mathrm{Spherical\ NFW}$')
	if triaxial: plt.hist(dist_tri,color='b',bins=nbins,histtype='step',normed=True,lw=3,weights=mass,label=r'$\mathrm{Triaxial\ NFW}$')
	if compare: plt.hist(dist_mt,color='g',bins=nbins,histtype='step',normed=True,lw=3,weights=mass,label=r'$\mathrm{Merger\ Tree}$')
	plt.title(r'$\mathrm{Host\ %i}$' %host_idx)
	plt.xlabel(r'$\mathrm{r/R_{200}}$')
	plt.ylabel(r'$\mathrm{mass\ weighted\ subhalo\ distribution}$')
	plt.grid()
	plt.xlim([0,5])
	plt.savefig('plots/subhalo_distribution_hist_%i_%s_%s_%.0e.png' %(host_idx,integrator,potential,dt))

	# plt.show()