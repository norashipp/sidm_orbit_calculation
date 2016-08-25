import numpy as np
import matplotlib.pyplot as plt
import cPickle
import sys

import astropy.units as u
from gala.units import UnitSystem
import gala.potential as gp
import gala.dynamics as gd

from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.utils.setup import *

dpi = 175
fontsize = 12
plt.rc('savefig', dpi=dpi)
plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.rc('xtick.major', pad=5)
plt.rc('xtick.minor', pad=5)
plt.rc('ytick.major', pad=5)
plt.rc('ytick.minor', pad=5)

host_idx = int(sys.argv[1])
potential = 'triaxial_NFW'
# potential = 'triaxial_NFW'
dt = 4e-3
subs = np.array(sys.argv[2:],dtype=int)

host = HostHalo(idx=host_idx, potential=potential)
host.update(host.cosmo.age(0)) # potential will not be evolving

gravity = GetGravitationalForce(host)

usys = UnitSystem(u.Mpc, u.Gyr, u.Msun, u.degree, u.Mpc/u.Gyr)

M = host.mass_function(host=host, a=0, b=host.R_s)
R = host.R_s
fg, _ = gravity.calculate_partial_force([host.R_s,0,0])
v = np.sqrt(fg*host.R_s)

if host.potential == 'spherical_NFW':
	pot = gp.SphericalNFWPotential(v_c = v*u.Mpc/u.Gyr, r_s=R*u.Mpc, units=usys)
if host.potential == 'triaxial_NFW':
	pot = gp.LeeSutoTriaxialNFWPotential(v_c =v*u.Mpc/u.Gyr, r_s = R*u.Mpc, a = 1, b = host.q, c = host.s, units=usys)

for sub_idx in subs:
	sub = host.subhalos[sub_idx]
	if sub:
		# host.update(sub.t0)
		time = host.cosmo.age(0) - sub.t0
		n_steps = time/dt

		w0 = gd.CartesianPhaseSpacePosition(pos=sub.position*u.Mpc,vel=sub.momentum*u.Mpc/u.Gyr)

		orbit = pot.integrate_orbit(w0, dt=dt, n_steps=n_steps)
		
		write = 1
		if write:
			f = open(HOMEDIR + '/sidm_orbit_calculation/src/output/gala_orbit_%i_%s_%.0e_%i.dat'%(host_idx,potential,dt,sub_idx),'wb')
			cPickle.dump(orbit,f)
			f.close()

		plot = 0
		if plot:
			print 'plotting...'
			fig = orbit.plot()
			plt.show()
	else:
		print sub_idx
