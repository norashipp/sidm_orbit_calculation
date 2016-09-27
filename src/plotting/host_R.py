import numpy as np
from numpy import savetxt

from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.utils.setup import *

R = []
potential = 'spherical_NFW'
for host_idx in range(51):
	host = HostHalo(host_idx,potential,subs=False,scale_density=False)
	host.update(host.cosmo.age(0))

	print 'Host %i' %host_idx
	print 'R = %.2e' %host.R

	R.append(host.R)

R = np.asarray(R)
np.savetxt(HOMEDIR+'sidm_orbit_calculation/src/plotting/host_R.txt',R)
