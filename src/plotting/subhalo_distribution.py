import cPickle
import numpy as np
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt

from sidm_orbit_calculation.src.utils.setup import *

def kde_sklearn(data, grid, bandwidth = 1.0, **kwargs):
    kde_skl = KernelDensity(bandwidth = bandwidth, **kwargs)
    kde_skl.fit(data[:, np.newaxis])
    log_pdf = kde_skl.score_samples(grid[:, np.newaxis]) # sklearn returns log(density)

    return np.exp(log_pdf)

f = open(HOMEDIR + 'sidm_orbit_calculation/src/final_distribution_40.dat','rb')
dat = cPickle.load(f)
f.close()

host_radius = 1.21100235447

sigma = np.std(dat)
n = len(dat)
h = 1.059*sigma*n**-(1/5)
print 'h = ', h

grid = np.arange(0,7,0.1)
kde = kde_sklearn(dat/host_radius, grid, h, kernel='tophat')

plt.figure()
plt.plot(grid,kde)
plt.title('Tophat KDE')
plt.xlabel('r/R_200')
plt.ylabel('subhalos')
plt.savefig('41_distribution_kde.png')

plt.figure()
plt.hist(dat,bins=30,histtype='step',normed=True)
plt.title('Histogram')
plt.xlabel('r/R_200')
plt.ylabel('subhalos')
plt.savefig('41_distribution_hist.png')
