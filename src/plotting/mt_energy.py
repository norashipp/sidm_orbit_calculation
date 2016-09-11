import numpy as np
import sys
import matplotlib.pyplot as plt

from colossus.cosmology import cosmology

from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.utils.setup import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.merger_tree.cluster import *
from sidm_orbit_calculation.src.calculation.get_gravitational_force import *

dpi = 175
fontsize = 15
plt.rc('savefig', dpi=dpi)
plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.rc('xtick.major', pad=5)
plt.rc('xtick.minor', pad=5)
plt.rc('ytick.major', pad=5)
plt.rc('ytick.minor', pad=5)

my_cosmo = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.045714, 'sigma8': 0.82, 'ns': 0.96}
cosmo = cosmology.setCosmology('my_cosmo', my_cosmo)

host_idx = int(sys.argv[1])
# sub_idx = int(sys.argv[2])
sub_idx_array = np.array(sys.argv[2:],dtype=int)
# sub_idx_array = np.arange(16)
integrator = 'leapfrog'
dt = 4e-3
potential = 'triaxial_NFW_BT'

subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)

host = HostHalo(host_idx,potential,subs=False)
# host.update(cosmo.age(0))
gravity = GetGravitationalForce(host)

for sub_idx in sub_idx_array:
	print 'Plotting subhalo %i' %sub_idx
	sub = subs[sub_idx]

	t_mt = cosmo.age(1/sub.a-1)
	
	h = 0.7
	x_mt = sub.rel_x/(h*(1/sub.a))
	y_mt = sub.rel_y/(h*(1/sub.a))
	z_mt = sub.rel_z/(h*(1/sub.a))

	vx_mt = sub.rel_vx*1000*m_to_Mpc/s_to_Gyr
	vy_mt = sub.rel_vy*1000*m_to_Mpc/s_to_Gyr
	vz_mt = sub.rel_vz*1000*m_to_Mpc/s_to_Gyr

	dist_mt = np.sqrt(x_mt**2 + y_mt**2 + z_mt**2)
	vt_mt = np.sqrt(vx_mt**2 + vy_mt**2 + vz_mt**2)

	energy_array = []
	for i in range(len(t_mt)):
		host.update(t_mt[i])
		KE = 0.5*host.M*vt_mt[i]**2
		PE = gravity.potential_function(x_mt[i],y_mt[i],z_mt[i])
		# print t_mt[i], KE, PE
		# print np.sign(KE), np.sign(PE)
		E = KE + PE
		energy_array.append(E)

	plt.figure()
	plt.plot(t_mt,energy_array)
	# plt.savefig()
plt.show()
