# import matplotlib as mpl
# mpl.use('Agg')

import numpy as np
import cPickle
import sys
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology

from sidm_orbit_calculation.src.utils.setup import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.merger_tree.cluster import *

my_cosmo = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.045714, 'sigma8': 0.82, 'ns': 0.96}
cosmo = cosmology.setCosmology('my_cosmo', my_cosmo)
h = 0.7

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
subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)

dthresh = 0.2
tmin = 12
tmax = 14

sub = subs[sub_idx_array[0]]
# t0 = cosmo.age(1/sub.a-1)
a0 = sub.a
x0 = sub.rel_x/(h*(1/sub.a))
y0 = sub.rel_y/(h*(1/sub.a))
z0 = sub.rel_z/(h*(1/sub.a))

a1 = a0[-6:]
# idx0 = np.where(np.all([t0>tmin,t0<tmax],axis=0))

fig, ax =  plt.subplots(1,3,figsize=(18,5))
ax[0].plot(x0[-6:],y0[-6:],'b-',lw=3)
ax[1].plot(y0[-6:],z0[-6:],'b-',lw=3)
ax[2].plot(z0[-6:],x0[-6:],'b-',lw=3)
ax[0].plot(x0[-6:],y0[-6:],'b.',ms=10)
ax[1].plot(y0[-6:],z0[-6:],'b.',ms=10)
ax[2].plot(z0[-6:],x0[-6:],'b.',ms=10)
ax[0].plot(x0[-4],y0[-4],'b^',ms=10)
ax[1].plot(y0[-4],z0[-4],'b^',ms=10)
ax[2].plot(z0[-4],x0[-4],'b^',ms=10)

ax[0].set_xlabel('x')
ax[0].set_ylabel('y')
ax[1].set_xlabel('y')
ax[1].set_ylabel('z')
ax[2].set_xlabel('z')
ax[2].set_ylabel('x')

count = 0
colors = np.array(['g','r','c','k','m','y'],dtype=str)
for sub_idx in range(len(subs)-1):
    if sub_idx == sub_idx_array[0]: continue
    sub = subs[sub_idx]
    if sub.v_max[-1] < 70: continue
    # t = cosmo.age(1/sub.a-1)
    a = sub.a
    x = sub.rel_x/(h*(1/sub.a))
    y = sub.rel_y/(h*(1/sub.a))
    z = sub.rel_z/(h*(1/sub.a))
    # d = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)

    xi = interp1d(a,x)
    yi = interp1d(a,y)
    zi = interp1d(a,z)

    # di = interp1d(a,d)

    x1 = xi(a1)
    y1 = yi(a1)
    z1 = zi(a1)
    d1 = np.sqrt((x1-x0[-6])**2 + (y1-y0[-6])**2 + (z1-z0[-6])**2)

    if np.any(d1 < dthresh):
        c = colors[count]
        # c = k
        count += 1
        ax[0].plot(x1,y1,'--',lw=3,ms=10,c=c)
        ax[1].plot(y1,z1,'--',lw=3,ms=10,c=c)
        ax[2].plot(z1,x1,'--',lw=3,ms=10,c=c)
        ax[0].plot(x1,y1,'.',lw=3,ms=10,c=c)
        ax[1].plot(y1,z1,'.',lw=3,ms=10,c=c)
        ax[2].plot(z1,x1,'.',lw=3,ms=10,c=c)
        ax[0].plot(xi(a1[-4]),yi(a1[-4]),'^',lw=3,ms=10,c=c)
        ax[1].plot(yi(a1[-4]),zi(a1[-4]),'^',lw=3,ms=10,c=c)
        ax[2].plot(zi(a1[-4]),xi(a1[-4]),'^',lw=3,ms=10,c=c)

print '%i subhalos within %.2f Mpc' %(count,dthresh)    
ax[1].set_title('%i subhalos within %.2f Mpc' %(count,dthresh))


    # print t.min(),t.max()
    # idx = np.where(np.all([t>tmin,t<tmax],axis=0))
    # x = x[idx]
    # y = y[idx]
    # z = z[idx]
    # d = np.sqrt((x-x0[idx0])**2 + (y-y0[idx0])**2 + (z-z0[idx0])**2)

    # di = interp1d(t,d)
    # idx = np.where(t>tmin)
    # if np.any(d < dthresh):
    #     ax[0].plot(x[idx],y[idx],'ko',lw=3,ms=10)
    #     ax[0].plot(y[idx],z[idx],'ko',lw=3,ms=10)
    #     ax[0].plot(z[idx],x[idx],'ko',lw=3,ms=10)

plt.show()
