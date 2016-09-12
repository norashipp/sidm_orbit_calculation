import numpy as np
import matplotlib.pyplot as plt

from sidm_orbit_calculation.src.merger_tree.cluster import *

dpi = 175
fontsize = 15
plt.rc('savefig', dpi=dpi)
# plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.rc('xtick.major', pad=5)
plt.rc('xtick.minor', pad=5)
plt.rc('ytick.major', pad=5)
plt.rc('ytick.minor', pad=5)

vmax = []
for i in range(51):
# for i in [12]:
	subs = SubHalos("merger_tree/subs/sub_%d.dat" %i)
	for j in range(len(subs)-1):
		vmax.append(subs[j].v_max[-1])
		print i, j, len(subs[j].v_max)

vmax = np.asarray(vmax)
print len(vmax)

print len(vmax)
bins = np.logspace(0,3)
hist, bins = np.histogram(vmax,bins)
cumsum = []
for i in range(len(hist)):
	cs = np.sum(hist[i:])
	cumsum.append(cs)

plt.figure()
plt.plot(bins[:-1],cumsum,lw=3)
plt.xlabel('vmax (km/s)')
plt.ylabel('N (>vmax)')
plt.grid()
plt.xscale('log')
plt.yscale('log')

plt.savefig('cumulative_vmax_dist.png')
plt.show()