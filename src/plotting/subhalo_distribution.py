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

'''
def kde_sklearn(data, grid, bandwidth = 1.0, **kwargs):
    kde_skl = KernelDensity(bandwidth = bandwidth, **kwargs)
    kde_skl.fit(data[:, np.newaxis])
    log_pdf = kde_skl.score_samples(grid[:, np.newaxis]) # sklearn returns log(density)

    return np.exp(log_pdf)
'''

for i in range(6):
	host_idx = i
	# host_idx = int(sys.argv[1])
	potential = 'spherical_NFW'

	dat = loadtxt(HOMEDIR + 'sidm_orbit_calculation/src/output/final_distances_%i.txt' %host_idx, unpack=True)

	# host_radius = 1.21100235447 # 40
	host = HostHalo(host_idx,potential,subs=False)
	host.update(host.cosmo.age(0))
	host_radius = host.R

	compare = 1
	if compare:
		subs = mt.SubHalos(HOMEDIR + 'sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat' % host_idx)
		merger_tree = []
		for i in range(len(subs)-1):
			sub = subs[i]
			dist = np.sqrt(sub.rel_x[-1]**2 + sub.rel_y[-1]**2 + sub.rel_z[-1]**2)
			merger_tree.append(dist/host_radius)

	sigma = np.std(dat)
	n = len(dat)
	h = 1.059*sigma*n**-(1/5)
	print 'h = ', h

	# grid = np.arange(0,20,0.1)
	# kde = kde_sklearn(dat/host_radius, grid, h, kernel='gaussian')

	'''
	plt.figure()
	plt.plot(grid,kde,lw=3)
	plt.title(r'$\mathrm{Host\ %i\ -\ Gaussian\ KDE}$' %host_idx)
	plt.xlabel(r'$\mathrm{r/R_{200}}$')
	plt.ylabel(r'$\mathrm{subhalos}$')
	plt.grid()
	# plt.savefig('%i_distribution_kde.png' %host_idx)
	'''

	# q75, q25 = np.percentile(dat, [75 ,25])
	# iqr = q75 - q25
	# h = 2*iqr*len(dat)**(-1/3)
	# nbins = (dat.max()-dat.min())/h
	# nbins = np.sqrt(len(dat))
	# print nbins
	nbins = 40

	plt.figure()
	plt.hist(dat,bins=nbins,histtype='step',normed=True,lw=3)
	if compare: plt.hist(merger_tree,bins=nbins,histtype='step',normed=True,lw=3)
	plt.title(r'$\mathrm{Host\ %i\ -\ Histogram}$' %host_idx)
	plt.xlabel(r'$\mathrm{r/R_{200}}$')
	plt.ylabel(r'$\mathrm{subhalos}$')
	plt.grid()
	plt.xlim([0,10])
	plt.savefig('plots/%i_distribution_hist.png' %host_idx)

	# plt.show()