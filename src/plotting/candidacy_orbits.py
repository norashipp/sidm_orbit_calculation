import matplotlib as mpl
mpl.use('Agg') # not interactive

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

my_cosmo = {'flat': True, 'H0': 70.0, 'Om0': 0.27, 'Ob0': 0.045714, 'sigma8': 0.82, 'ns': 0.96}
cosmo = cosmology.setCosmology('my_cosmo', my_cosmo)

dpi = 175
fontsize = 20
plt.rc('savefig', dpi=dpi)
plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.rc('xtick.major', pad=5)
plt.rc('xtick.minor', pad=5)
plt.rc('ytick.major', pad=5)
plt.rc('ytick.minor', pad=5)

host_idx = int(sys.argv[1])
sub_idx_array = np.array(sys.argv[2:],dtype=int)
integrator = 'leapfrog'
potential = 'spherical_NFW'
dt = 4e-3

host = HostHalo(host_idx,potential)



subs = SubHalos(HOMEDIR + "sidm_orbit_calculation/src/merger_tree/subs/sub_%d.dat" % host_idx)

c_mt = 'g'
c = 'r'
c_d = 'b'

ls_mt = '--'
ls = '-'
ls_d = '-'

for sub_idx in sub_idx_array:
    if host.subhalos[sub_idx]:
            fig1, ax1 = plt.subplots(1,1,figsize=(8,8))
            fig2, ax2 = plt.subplots(1,1,figsize=(8,8))
            fig3, ax3 = plt.subplots(1,1,figsize=(8,8))
            # fig.tight_layout(pad=5.0, w_pad=2.0, h_pad=2.0)
            # fig.set_tight_layout(True)
            rmax = 1

            print 'Plotting subhalo %i' %sub_idx
            sub = subs[sub_idx]

            ### merger tree ###
            merger = 1
            if merger:
                    h = 0.7
                    
                    t_mt = cosmo.age(1/sub.a-1)
                    idx = t_mt > host.subhalos[sub_idx].t0
                    t_mt = t_mt[idx]

                    x_mt = sub.rel_x/(h*(1/sub.a))
                    y_mt = sub.rel_y/(h*(1/sub.a))
                    z_mt = sub.rel_z/(h*(1/sub.a))

                    vx_mt = sub.rel_vx # *1000*m_to_Mpc/s_to_Gyr
                    vy_mt = sub.rel_vy # *1000*m_to_Mpc/s_to_Gyr
                    vz_mt = sub.rel_vz # *1000*m_to_Mpc/s_to_Gyr

                    x_mt = x_mt[idx]
                    y_mt = y_mt[idx]
                    z_mt = z_mt[idx]

                    vx_mt = vx_mt[idx]
                    vy_mt = vy_mt[idx]
                    vz_mt = vz_mt[idx]

                    radii_mt = []
                    for time in t_mt:
                        if time <= host.cosmo.age(0):
                            # print time
                            host.update(time)
                            radii_mt.append(host.R)
                    radii_mt = np.asarray(radii_mt)
                    
                    dist_mt = np.sqrt(x_mt**2 + y_mt**2 + z_mt**2)/radii_mt
                    vt_mt = np.sqrt(vx_mt**2 + vy_mt**2 + vz_mt**2)

                    x_mt/=host.R
                    y_mt/=host.R
                    z_mt/=host.R

                    ax1.plot(x_mt, y_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    ax2.plot(t_mt, dist_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    ax3.plot(0,0, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')

            ### no drag ###

            infile = HOMEDIR+'data/candidacy/sigma0/%i_%s_%.0e_0.00_%i.dat' %(host_idx,potential,dt,sub_idx)
            # infile = HOMEDIR+'sidm_orbit_calculation/src/output/%i_leapfrog_%s_%.0e_%i.dat' %(host_idx,potential,dt,sub_idx)
            
            f = open(infile,'rb')
            data = cPickle.load(f)
            f.close()
            t, positions, velocities = data

            x = positions[:,0]
            y = positions[:,1]
            z = positions[:,2]

            vx = velocities[:,0]/(1000*m_to_Mpc/s_to_Gyr)
            vy = velocities[:,1]/(1000*m_to_Mpc/s_to_Gyr)
            vz = velocities[:,2]/(1000*m_to_Mpc/s_to_Gyr)

            print len(t), len(x)

            radii = []
            for time in t:
                if time > host.cosmo.age(0): time = host.cosmo.age(0)
                # print time
                host.update(time)
                radii.append(host.R)
            radii = np.asarray(radii)
            host.update(host.cosmo.age(0))
            
            dist = np.sqrt(x**2 + y**2 + z**2)/radii
            vt = np.sqrt(vx**2 + vy**2 + vz**2)

            x/=host.R
            y/=host.R
            z/=host.R

            sigma = 0.
            ax1.plot(x, y, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            ax2.plot(t,dist, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            ax3.plot(0,0, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            
            # sigs = [3,6,9,12,15,18]
            sigs = [3,9,15,21]
            for sigma in sigs:
                    infile_drag = HOMEDIR+'data/candidacy/sigma%i/%i_%s_%.0e_%.2f_%i.dat' %(sigma,host_idx,potential,dt,sigma,sub_idx)
                    # infile_drag = HOMEDIR+'sidm_orbit_calculation/src/output/sigma%i/%i_dissipative_%s_%.0e_%i.dat' %(sigma,host_idx,potential,dt,sub_idx)
                    f = open(infile_drag,'rb')
                    data = cPickle.load(f)
                    f.close()
                    _, positions_drag, velocities_drag = data
                    x_d = positions_drag[:,0]
                    y_d = positions_drag[:,1]
                    z_d = positions_drag[:,2]

                    dist_d = np.sqrt(x_d**2 + y_d**2 + z_d**2)/radii

                    vx_d = velocities_drag[:,0]/(1000*m_to_Mpc/s_to_Gyr)
                    vy_d = velocities_drag[:,1]/(1000*m_to_Mpc/s_to_Gyr)
                    vz_d = velocities_drag[:,2]/(1000*m_to_Mpc/s_to_Gyr)

                    vt_d = np.sqrt(vx_d**2 + vy_d**2 + vz_d**2)

                    x_d/=host.R
                    y_d/=host.R
                    z_d/=host.R
                    
                    ax1.plot(x_d, y_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    ax2.plot(t, dist_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    ax3.plot(0, 0, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    
            triaxial = 0
            if triaxial:
                    integrator = 'leapfrog'
                    potential = 'triaxial_NFW_BT'
                    infile_triaxial = HOMEDIR+'sidm_orbit_calculation/src/output/%i_%s_%s_%.0e_%i_major_axis.dat' %(host_idx,integrator,potential,dt,sub_idx)
                    # infile_triaxial = HOMEDIR+'sidm_orbit_calculation/src/output/%i_%s_spherical_NFW_%.0e_%i_major_axis.dat' %(host_idx,integrator,dt,sub_idx)
                    f = open(infile_triaxial,'rb')
                    data = cPickle.load(f)
                    f.close()
                    _, positions_tri, velocities_tri = data
                    x_tri = positions_tri[:,0]/host.R
                    y_tri = positions_tri[:,1]/host.R
                    z_tri = positions_tri[:,2]/host.R

                    dist_tri = np.sqrt(x_tri**2 + y_tri**2 + z_tri**2)

                    vx_tri = velocities_tri[:,0]/(1000*m_to_Mpc/s_to_Gyr)
                    vy_tri = velocities_tri[:,1]/(1000*m_to_Mpc/s_to_Gyr)
                    vz_tri = velocities_tri[:,2]/(1000*m_to_Mpc/s_to_Gyr)

                    vt_tri = np.sqrt(vx_tri**2 + vy_tri**2 + vz_tri**2)

                    ax1.plot(x_tri, y_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
                    ax2.plot(t, dist_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
                    ax3.plot(0,0, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')

            ####### PLOTTING #######

            ####### ORBITS #######

            ax1.plot(0,0,'k*',markersize=12) # , label=r'$\mathrm{Host\ Center}$')
            ax1.plot(x[0],y[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # , label=r'$\mathrm{Orbit\ Start}$')
            ax1.set_xlabel(r'$x\ (R_{\rm 200m})$')
            ax1.set_ylabel(r'$y\ (R_{\rm 200m})$')
            
            ax2.plot(t[0],dist[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
            ax2.set_xlabel(r'$t\ {\rm ( Gyr)}$')
            ax2.set_ylabel(r'$r\ (R_{\rm 200m})$')

            ax3.legend(loc='center',fontsize=30)
            ax3.get_xaxis().set_visible(False)
            ax3.get_yaxis().set_visible(False)
            ax3.spines['bottom'].set_color('white')
            ax3.spines['top'].set_color('white') 
            ax3.spines['right'].set_color('white')
            ax3.spines['left'].set_color('white')
            
            fig1.savefig(HOMEDIR + '/sidm_orbit_calculation/src/plots/%i_%.0e_%i_x_y.png'%(host_idx,dt,sub_idx))
            fig2.savefig(HOMEDIR + '/sidm_orbit_calculation/src/plots/%i_%.0e_%i_radius.png'%(host_idx,dt,sub_idx))
            fig3.savefig(HOMEDIR + '/sidm_orbit_calculation/src/plots/%i_%.0e_%i_.png'%(host_idx,dt,sub_idx))

# plt.show()
