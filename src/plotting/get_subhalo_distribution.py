from __future__ import division
import pylab
import numpy as np
import cPickle
import sys
import glob

import matplotlib.pyplot as plt

from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.utils.setup import *
from sidm_orbit_calculation.src.halos.host_halo import *

host_idx = int(sys.argv[1])
integrator = 'leapfrog'
potential = 'spherical_NFW'

infiles = glob.glob(HOMEDIR+'sidm_orbit_calculation/src/output/%i_%s_%s_*.dat' %(host_idx,integrator,potential))
infiles.sort()

# final_positions = np.array([])
final_dists = []

# host = HostHalo(host_idx,potential)
# host.update(host.cosmo.age(0))
# print host.R
host_radius = 1.21100235447
i = 0
for infile in infiles:
	print i
	f = open(infile,'rb')
	data = cPickle.load(f)
	f.close()
	# times,positions,momenta,gravity,drag,density,energy,host_idx = data
	times,positions,momenta,gravity,drag,density,energy,host_idx,potential,host_radius = data

	fp = positions[-1]
	dist = np.sqrt(fp[0]**2 + fp[1]**2  + fp[2]**2)
	final_dists.append(dist)
	i+=1
	
final_dists = np.asarray(final_dists)

f = open('final_distribution_40.dat','wb')
cPickle.dump(final_dists,f)
f.close()

plt.hist(final_dists/host_radius,bins=50)
plt.show()
