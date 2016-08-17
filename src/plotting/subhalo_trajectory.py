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
times,positions,momenta,gravity,drag,density,energy,host_idx,potential,host_radius = data

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
lim = 2*host_radius

# mt_x = mt_x[mt_dist < lim]
# mt_y = mt_y[mt_dist < lim]
# mt_z = mt_z[mt_dist < lim]

# print x.min(), x.max()
# print mt_x.min(), mt_x.max()

print mt_x[0], mt_y[0]
print x[0], y[0]

plt.figure()
plt.plot(x, y, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
plt.plot(mt_x, mt_y, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
plt.plot(0,0,'^k',markersize=12,label=r'$\mathrm{Host\ Center}$')
plt.plot(x[0],y[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
# plt.plot(mt_x[0],mt_y[0],'g*',markersize=10)
# plt.plot(1.03121194,-0.90783956,'r*',markersize=10) # 40, 100
# plt.plot(-1.29679797, -0.45553861,'r*',markersize=10) # 40, 0
# plt.plot([-host_radius,host_radius],[-0.45553861,-0.45553861],'r--',lw=3)
plt.xlabel(r'$\mathrm{x\ (Mpc)}$')
plt.ylabel(r'$\mathrm{y\ (Mpc)}$')
plt.title(r'$\mathrm{Host\ %i,\ Subhalo\ %i}$' %(host_idx, sub_idx))
plt.legend(loc='lower left',fontsize=15)
plt.xlim([-2,2])
plt.ylim([-2,2])
plt.grid()

plt.show()
