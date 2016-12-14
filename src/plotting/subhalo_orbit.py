import matplotlib as mpl
# mpl.use('Agg')

import numpy as np
import cPickle
import sys
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.cm import ocean

from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.utils.setup import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.merger_tree.cluster import *

my_cosmo = {'flat': True, 'H0': 70.0, 'Om0': 0.27,
            'Ob0': 0.045714, 'sigma8': 0.82, 'ns': 0.96}
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
# sub_idx = int(sys.argv[2])
sub_idx_array = np.array(sys.argv[2:], dtype=int)
# sub_idx_array = np.arange(21)
integrator = 'leapfrog'
dt = 4e-3
potential = 'spherical_NFW'
host = HostHalo(host_idx, potential)
host.update(host.cosmo.age(0))

subs = SubHalos(DATADIR + '/subs/sub_%d.dat' % host_idx)

for sub_idx in sub_idx_array:
    if host.subhalos[sub_idx]:
        print 'Plotting subhalo %i' % sub_idx
        potential = 'spherical_NFW'
        infile = HOMEDIR + \
            'scratch-midway/%i_%s_%.0e_0.00_%i.dat' % (
                host_idx, potential, dt, sub_idx)
        sub = subs[sub_idx]

        f = open(infile, 'rb')
        data = cPickle.load(f)
        f.close()
        t, positions, velocities = data

        x = positions[:, 0]
        y = positions[:, 1]
        z = positions[:, 2]

        vx = velocities[:, 0]
        vy = velocities[:, 1]
        vz = velocities[:, 2]

        print 'initial velocities = ', vx[0], vy[0], vz[0]

        vx /= (1000 * m_to_Mpc / s_to_Gyr)
        vy /= (1000 * m_to_Mpc / s_to_Gyr)
        vz /= (1000 * m_to_Mpc / s_to_Gyr)

        radii = []
        for time in t:
            if time < host.cosmo.age(0):
                # print time
                host.update(time)
            radii.append(host.R)

        radii = np.asarray(radii)

        dist = np.sqrt(x**2 + y**2 + z**2)
        dist /= radii
        vt = np.sqrt(vx**2 + vy**2 + vz**2)

        host.update(host.cosmo.age(0))
        x /= host.R
        y /= host.R
        z /= host.R

        ### merger tree ###

        t_mt = cosmo.age(1 / sub.a - 1)
        # idx = t_mt >= host.subhalos[sub_idx].t0
        # if idx[0] > 0: idx = np.insert(idx,0,idx[0]-1)
        idx = np.arange(len(t_mt))
        idx[np.where(idx == True)[0][0] - 1] = True
        t_mt = t_mt[idx]

        h = 0.7
        x_mt = sub.rel_x / (h * (1 / sub.a))
        y_mt = sub.rel_y / (h * (1 / sub.a))
        z_mt = sub.rel_z / (h * (1 / sub.a))

        vx_mt = sub.rel_vx  # *1000*m_to_Mpc/s_to_Gyr
        vy_mt = sub.rel_vy  # *1000*m_to_Mpc/s_to_Gyr
        vz_mt = sub.rel_vz  # *1000*m_to_Mpc/s_to_Gyr

        x_mt = x_mt[idx]
        y_mt = y_mt[idx]
        z_mt = z_mt[idx]

        vx_mt = vx_mt[idx]
        vy_mt = vy_mt[idx]
        vz_mt = vz_mt[idx]

        dist_mt = np.sqrt(x_mt**2 + y_mt**2 + z_mt**2)
        vt_mt = np.sqrt(vx_mt**2 + vy_mt**2 + vz_mt**2)

        radii = []
        for time in t_mt:
            if time < host.cosmo.age(0):
                # print time
                host.update(time)
            radii.append(host.R)

        radii = np.asarray(radii)
        dist_mt /= radii
        host.update(host.cosmo.age(0))

        x_mt /= host.R
        y_mt /= host.R
        z_mt /= host.R

        fig, ax = plt.subplots(3, 3, figsize=(20, 20))
        fig.tight_layout(pad=5.0, w_pad=2.0, h_pad=2.0)

        ####### ORBITS #######

        c = 'slateblue'
        c_mt = 'forestgreen'

        ls = '-'
        ls_mt = '--'

        # rmax = 1

        ax[0][0].plot(x_mt, y_mt, color=c_mt, ls=ls_mt,
                      lw=3, label=r'$\mathrm{Merger\ Tree}$')
        ax[0][0].plot(x, y, color=c, ls=ls, lw=3,
                      label=r'$\mathrm{Spherical\ NFW}$')
        # , label=r'$\mathrm{Host\ Center}$')
        ax[0][0].plot(0, 0, 'k*', markersize=12)
        # , label=r'$\mathrm{Orbit\ Start}$')
        ax[0][0].plot(x[0], y[0], '*', color=c, ls=ls,
                      markeredgecolor=None, markersize=15)
        ax[0][0].set_xlabel(r'$\mathrm{x\ (R_{200m})}$')
        ax[0][0].set_ylabel(r'$\mathrm{y\ (R_{200m})}$')

        ax[0][1].plot(y_mt, z_mt, color=c_mt, ls=ls_mt,
                      lw=3, label=r'$\mathrm{Merger\ Tree}$')
        ax[0][1].plot(y, z, color=c, ls=ls, lw=3,
                      label=r'$\mathrm{Spherical\ NFW}$')
        # , label=r'$\mathrm{Host\ Center}$')
        ax[0][1].plot(0, 0, 'k*', markersize=12)
        # ,label=r'$\mathrm{Orbit\ Start}$')
        ax[0][1].plot(y[0], z[0], '*', color=c, ls=ls,
                      markeredgecolor=None, markersize=15)
        ax[0][1].set_xlabel(r'$\mathrm{y\ (R_{200m})}$')
        ax[0][1].set_ylabel(r'$\mathrm{z\ (R_{200m})}$')

        ax[0][2].plot(z_mt, x_mt, color=c_mt, ls=ls_mt,
                      lw=3, label=r'$\mathrm{Merger\ Tree}$')
        ax[0][2].plot(z, x, color=c, ls=ls, lw=3,
                      label=r'$\mathrm{Spherical\ NFW}$')
        # , label=r'$\mathrm{Host\ Center}$')
        ax[0][2].plot(0, 0, 'k*', markersize=12)
        # , label=r'$\mathrm{Orbit\ Start}$')
        ax[0][2].plot(z[0], x[0], '*', color=c, ls=ls,
                      markeredgecolor=None, markersize=15)
        ax[0][2].set_xlabel(r'$\mathrm{z\ (R_{200m})}$')
        ax[0][2].set_ylabel(r'$\mathrm{x\ (R_{200m})}$')

        ####### VELOCITIES #######

        ax[1][0].plot(t_mt, vx_mt, color=c_mt, ls=ls_mt,
                      lw=3, label=r'$\mathrm{Merger\ Tree}$')
        ax[1][0].plot(t, vx, color=c, ls=ls, lw=3,
                      label=r'$\mathrm{Spherical\ NFW}$')
        # ,label=r'$\mathrm{Orbit\ Start}$')
        ax[1][0].plot(t[0], vx[0], '*', color=c, ls=ls,
                      markeredgecolor=None, markersize=15)
        ax[1][0].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
        ax[1][0].set_ylabel(r'$\mathrm{vx\ (km/s)}$')

        ax[1][1].plot(t_mt, vy_mt, color=c_mt, ls=ls_mt,
                      lw=3, label=r'$\mathrm{Merger\ Tree}$')
        ax[1][1].plot(t, vy, color=c, ls=ls, lw=3,
                      label=r'$\mathrm{Spherical\ NFW}$')
        # ,label=r'$\mathrm{Orbit\ Start}$')
        ax[1][1].plot(t[0], vy[0], '*', color=c, ls=ls,
                      markeredgecolor=None, markersize=15)
        ax[1][1].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
        ax[1][1].set_ylabel(r'$\mathrm{vy\ (km/s)}$')

        ax[1][2].plot(t_mt, vz_mt, color=c_mt, ls=ls_mt,
                      lw=3, label=r'$\mathrm{Merger\ Tree}$')
        ax[1][2].plot(t, vz, color=c, ls=ls, lw=3,
                      label=r'$\mathrm{Spherical\ NFW}$')
        # ,label=r'$\mathrm{Orbit\ Start}$')
        ax[1][2].plot(t[0], vz[0], '*', color=c, ls=ls,
                      markeredgecolor=None, markersize=15)
        ax[1][2].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
        ax[1][2].set_ylabel(r'$\mathrm{vz\ (km/s)}$')

        ####### VELOCITY MAGNITUDE #######

        ax[2][0].plot(t_mt, vt_mt, color=c_mt, ls=ls_mt,
                      lw=3, label=r'$\mathrm{Merger\ Tree}$')
        ax[2][0].plot(t, vt, color=c, ls=ls, lw=3,
                      label=r'$\mathrm{Spherical\ NFW}$')
        # ,label=r'$\mathrm{Orbit\ Start}$')
        ax[2][0].plot(t[0], vt[0], '*', color=c, ls=ls,
                      markeredgecolor=None, markersize=15)
        ax[2][0].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
        ax[2][0].set_ylabel(r'$\mathrm{v\ (km/s)}$')

        ####### RADIUS #######

        ax[2][1].plot(t_mt, dist_mt, color=c_mt, ls=ls_mt,
                      lw=3, label=r'$\mathrm{Merger\ Tree}$')
        ax[2][1].plot(t, dist, color=c, ls=ls, lw=3,
                      label=r'$\mathrm{Spherical\ NFW}$')
        # ,label=r'$\mathrm{Orbit\ Start}$')
        ax[2][1].plot(t[0], dist[0], '*', color=c, ls=ls,
                      markeredgecolor=None, markersize=15)
        ax[2][1].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
        ax[2][1].set_ylabel(r'$\mathrm{r\ (R_{200m})}$')
        ax[2][1].set_yscale('log')

        ax[2][2].plot(0, 0, color=c_mt, ls=ls_mt, lw=3,
                      label=r'$\mathrm{Merger\ Tree}$')
        ax[2][2].plot(0, 0, color=c, ls=ls, lw=3,
                      label=r'$\mathrm{Spherical\ NFW}$')
        ax[2][2].legend(loc='center', fontsize=30)
        ax[2][2].get_xaxis().set_visible(False)
        ax[2][2].get_yaxis().set_visible(False)
        ax[2][2].spines['bottom'].set_color('white')
        ax[2][2].spines['top'].set_color('white')
        ax[2][2].spines['right'].set_color('white')
        ax[2][2].spines['left'].set_color('white')

        ax[0][1].set_title(r'$\mathrm{Host\ %i,\ Subhalo\ %i}$' % (
            host_idx, sub_idx), fontsize=30)

        plt.savefig(
            HOMEDIR + '/sidm_orbit_calculation/src/plots/subhalo_%i_%i.png' % (host_idx, sub_idx))

# plt.show()
