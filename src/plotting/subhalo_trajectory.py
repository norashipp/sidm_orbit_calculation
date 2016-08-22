import numpy as np
import cPickle
import sys
import matplotlib.pyplot as plt

from sidm_orbit_calculation.src.plotting.make_plots import *
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

host_idx = int(sys.argv[1])
sub_idx = int(sys.argv[2])
integrator = 'leapfrog'
potential = 'spherical_NFW'
dt = 4e-3

infile = HOMEDIR+'sidm_orbit_calculation/src/output/%i_%s_%s_%.0e_%i.dat' %(host_idx,integrator,potential,dt,sub_idx)

# host_radius = 1.21100235447 # host 40

f = open(infile,'rb')
data = cPickle.load(f)
f.close()
# times, positions, velocities = data
positions = data[1]

subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)

x = positions[:,0]
y = positions[:,1]
z = positions[:,2]
	
dist = np.sqrt(x**2 + y**2 + z**2)

sub = subs[sub_idx]
mt_x = sub.rel_x
mt_y = sub.rel_y
mt_z = sub.rel_z

mt_dist = np.sqrt(mt_x**2 + mt_y**2 + mt_z**2)

# GALA COMPARISON
gala = 1
if gala:
	f = open(HOMEDIR + 'sidm_orbit_calculation/src/output/gala_orbit_%i_%i.dat'%(host_idx,sub_idx))
	orbit = cPickle.load(f)
	f.close()
	x_g = orbit.w()[0]
	y_g = orbit.w()[1]
	z_g = orbit.w()[2]

# DRAG COMPARISON
compare = 0
if compare:
	integrator = 'dissipative'
	infile_drag = HOMEDIR+'sidm_orbit_calculation/src/output/%i_%s_%s_%.0e_%i.dat' %(host_idx,integrator,potential,dt,sub_idx)
	f = open(infile_drag,'rb')
	data = cPickle.load(f)
	f.close()
	# times,positions,momenta,gravity,drag,density,energy,host_idx,host_radius = data
	# times,positions,momenta,gravity,drag,density,energy,host_idx,potential,host_radius = data
	_, positions_drag = data
	x_d = positions_drag[:,0]
	y_d = positions_drag[:,1]
	z_d = positions_drag[:,2]

plt.figure()
plt.plot(x, y, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
plt.plot(mt_x, mt_y, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
plt.plot(0,0,'^k',markersize=12,label=r'$\mathrm{Host\ Center}$')
plt.plot(x[0],y[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if compare: plt.plot(x_d, y_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
if gala: plt.plot(x_g, y_g, '--', color='violet', lw=3, label=r'$\mathrm{Gala}$')
plt.xlabel(r'$\mathrm{x\ (Mpc)}$')
plt.ylabel(r'$\mathrm{y\ (Mpc)}$')
plt.title(r'$\mathrm{Host\ %i,\ Subhalo\ %i}$' %(host_idx, sub_idx))
plt.legend(loc='lower left',fontsize=15)
plt.xlim([-2,2])
plt.ylim([-2,2])
plt.grid()
# plt.savefig('%i_%i_trajectory_y_x.png' %(host_idx,sub_idx))

plt.figure()
plt.plot(y, z, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
plt.plot(mt_y, mt_z, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
plt.plot(0,0,'^k',markersize=12,label=r'$\mathrm{Host\ Center}$')
plt.plot(y[0],z[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if compare: plt.plot(y_d, z_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
if gala: plt.plot(y_g, z_g, '--', color='violet', lw=3, label=r'$\mathrm{Gala}$')
plt.xlabel(r'$\mathrm{y\ (Mpc)}$')
plt.ylabel(r'$\mathrm{z\ (Mpc)}$')
plt.title(r'$\mathrm{Host\ %i,\ Subhalo\ %i}$' %(host_idx, sub_idx))
plt.legend(loc='lower left',fontsize=15)
plt.xlim([-2,2])
plt.ylim([-2,2])
plt.grid()
# plt.savefig('%i_%i_trajectory_z_y.png' %(host_idx,sub_idx))


plt.figure()
plt.plot(z, x, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
plt.plot(mt_z, mt_x, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
plt.plot(0,0,'^k',markersize=12,label=r'$\mathrm{Host\ Center}$')
plt.plot(z[0],x[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if compare: plt.plot(z_d, x_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
if gala: plt.plot(z_g, x_g, '--', color='violet', lw=3, label=r'$\mathrm{Gala}$')
plt.xlabel(r'$\mathrm{z\ (Mpc)}$')
plt.ylabel(r'$\mathrm{x\ (Mpc)}$')
plt.title(r'$\mathrm{Host\ %i,\ Subhalo\ %i}$' %(host_idx, sub_idx))
plt.legend(loc='lower left',fontsize=15)
plt.xlim([-2,2])
plt.ylim([-2,2])
plt.grid()
# plt.savefig('%i_%i_trajectory_x_z.png' %(host_idx,sub_idx))

plt.show()
