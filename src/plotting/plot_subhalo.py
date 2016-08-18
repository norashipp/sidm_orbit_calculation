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
# times,positions,momenta,gravity,drag,density,energy,host_idx,host_radius = data
# times,positions,momenta,gravity,drag,density,energy,host_idx,potential,host_radius = data
times, positions, velocities = data

subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)

t = times[:-1]

x = positions[:,0]
y = positions[:,1]
z = positions[:,2]

vx = velocities[:,0]
vy = velocities[:,1]
vz = velocities[:,2]
	
dist = np.sqrt(x**2 + y**2 + z**2)

sub = subs[sub_idx]
mt_x = sub.rel_x[-len(t):]
mt_y = sub.rel_y[-len(t):]
mt_z = sub.rel_z[-len(t):]

mt_vx = sub.rel_vx[-len(t):]
mt_vy = sub.rel_vy[-len(t):]
mt_vz = sub.rel_vz[-len(t):]

mt_dist = np.sqrt(mt_x**2 + mt_y**2 + mt_z**2)

print mt_x[0], mt_y[0]
print x[0], y[0]

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


fig, ax = plt.subplots(3,3,sharey=True)
plt.title(r'$\mathrm{Host\ %i,\ Subhalo\ %i}$' %(host_idx, sub_idx))

ax[0].plot(x, y, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
ax[0].plot(mt_x, mt_y, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[0].plot(0,0,'^k',markersize=12,label=r'$\mathrm{Host\ Center}$')
ax[0].plot(x[0],y[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if compare: ax[0].plot(x_d, y_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
ax[0].set_xlabel(r'$\mathrm{x\ (Mpc)}$')
ax[0].set_ylabel(r'$\mathrm{y\ (Mpc)}$')
ax[0].legend(loc='lower left',fontsize=15)
ax[0].xlim([-2,2])
ax[0].ylim([-2,2])

ax[1].plot(y, z, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
ax[1].plot(mt_y, mt_z, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[1].plot(0,0,'^k',markersize=12,label=r'$\mathrm{Host\ Center}$')
ax[1].plot(y[0],z[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if compare: ax[1].plot(y_d, z_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
ax[1].xlabel(r'$\mathrm{y\ (Mpc)}$')
ax[1].ylabel(r'$\mathrm{z\ (Mpc)}$')
ax[1].legend(loc='lower left',fontsize=15)
ax[1].set_xlim([-2,2])
ax[1].set_ylim([-2,2])

ax[2].plot(z, x, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
ax[2].plot(mt_z, mt_x, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[2].plot(0,0,'^k',markersize=12,label=r'$\mathrm{Host\ Center}$')
ax[2].plot(z[0],x[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if compare: ax[2].plot(z_d, x_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
ax[2].set_xlabel(r'$\mathrm{z\ (Mpc)}$')
ax[2].set_ylabel(r'$\mathrm{x\ (Mpc)}$')
ax[2].legend(loc='lower left',fontsize=15)
ax[2].xlim([-2,2])
ax[2].ylim([-2,2])

plt.show()
