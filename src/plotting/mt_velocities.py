import matplotlib as mpl
mpl.use('Agg')

import numpy as np
import sys
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from scipy.interpolate import UnivariateSpline

from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.utils.setup import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.merger_tree.cluster import *

my_cosmo = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.045714, 'sigma8': 0.82, 'ns': 0.96}
cosmo = cosmology.setCosmology('my_cosmo', my_cosmo)

dpi = 175
fontsize = 15
plt.rc('savefig', dpi=dpi)
plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.rc('xtick.major', pad=5)
plt.rc('xtick.minor', pad=5)
plt.rc('ytick.major', pad=5)
plt.rc('ytick.minor', pad=5)

host_idx = int(sys.argv[1])
sub_idx_array = np.array(sys.argv[2:],dtype=int)
# sub_idx_array = np.arange(16)

subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)

for sub_idx in sub_idx_array:
	print 'Subhalo %i' %sub_idx
	sub = subs[sub_idx]

	t = cosmo.age(1/sub.a-1)
	
	h = 0.7
	x = sub.rel_x/(h*(1/sub.a))
	y = sub.rel_y/(h*(1/sub.a))
	z = sub.rel_z/(h*(1/sub.a))

	vx = sub.rel_vx # *1000*m_to_Mpc/s_to_Gyr
	vy = sub.rel_vy # *1000*m_to_Mpc/s_to_Gyr
	vz = sub.rel_vz # *1000*m_to_Mpc/s_to_Gyr

	dist = np.sqrt(x**2 + y**2 + z**2)
	vt = np.sqrt(vx**2 + vy**2 + vz**2)

	'''
	cvx = x[0]
	cvy = y[0]
	cvz = z[0]
	cumvx = [cvx]
	cumvy = [cvy]
	cumvz = [cvz]
	print cvx
	for i in range(len(t)-1):
		dt = t[i+1]-t[i]
		cvx += vx[i]*dt*1000*m_to_Mpc/s_to_Gyr
		cvy += vy[i]*dt*1000*m_to_Mpc/s_to_Gyr
		cvz += vz[i]*dt*1000*m_to_Mpc/s_to_Gyr
		print dt, cvx, vx[i]
		
		cumvx.append(cvx)
		cumvy.append(cvy)
		cumvz.append(cvz)
	'''

	print 'merger tree initial conditions'
	print 'x0 = ', x[0], y[0], z[0]
	print 'vz = ', vx[:5]
	print 'vy = ', vy[:5]
	print 'vz = ', vz[:5]

	s = 0.05
	xsp = UnivariateSpline(t,x,s=s)
	ysp = UnivariateSpline(t,y,s=s)
	zsp = UnivariateSpline(t,z,s=s)

	vxsp = xsp.derivative
	vysp = ysp.derivative
	vzsp = zsp.derivative

	plt.figure(figsize=(5,5))
	plt.plot(t,x,'.',markersize=2)
	plt.plot(t,xsp(t))
	plt.savefig('merger_tree_x_y_12_14.png')

	plt.figure(figsize=(5,5))
	plt.plot(t,y,'.',markersize=2)
	plt.plot(t,ysp(t))
	plt.savefig('merger_tree_y_z_12_14.png')

	plt.figure(figsize=(5,5))
	plt.plot(t,z,'.',markersize=2)
	plt.plot(t,zsp(t))
	plt.savefig('merger_tree_z_x_12_14.png')

	'''
	plt.figure(figsize=(5,5))
	plt.plot(cumvx,cumvy,lw=3)
	plt.plot(cumvx[0]/host.R,cumvy[0]/host.R,'*',markersize=10)
	plt.savefig('cumulative_velocity_x_y_12_14.png')

	plt.figure(figsize=(5,5))
	plt.plot(cumvy,cumvz,lw=3)
	plt.plot(cumvy[0]/host.R,cumvz[0]/host.R,'*',markersize=10)
	plt.savefig('cumulative_velocity_y_z_12_14.png')

	plt.figure(figsize=(5,5))
	plt.plot(cumvz,cumvx,lw=3)
	plt.plot(cumvz[0]/host.R,cumvx[0]/host.R,'*',markersize=10)
	plt.savefig('cumulative_velocity_z_x_12_14.png')
	'''
