from __future__ import division
import pylab
import numpy as np
import cPickle
import sys
import matplotlib.pyplot as plt

from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.utils.setup import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.merger_tree.cluster import *

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
plt.plot(x, y, 'c', lw=2, label='orbit calculation')
plt.plot(mt_x, mt_y, 'g', lw=2, label='merger tree')
plt.plot(0,0,'^k')
plt.plot(x[0],y[0],'c*',markersize=12)
plt.plot(mt_x[0],mt_y[0],'g*',markersize=10)
# plt.plot(1.03121194,-0.90783956,'r*',markersize=10) # 40, 100
# plt.plot(-1.29679797, -0.45553861,'r*',markersize=10) # 40, 0
# plt.plot([-host_radius,host_radius],[-0.45553861,-0.45553861],'r--',lw=3)
plt.xlabel('x (Mpc)')
plt.ylabel('y (Mpc)')

plt.show()
