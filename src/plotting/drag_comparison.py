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
            fig, ax = plt.subplots(3,3,figsize=(20,20))
            # fig1, ax1 = plt.subplots(1,3,figsize=(40,15))
            fig.tight_layout(pad=5.0, w_pad=2.0, h_pad=2.0)
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

                    ax[0][0].plot(x_mt, y_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    ax[0][1].plot(y_mt, z_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    ax[0][2].plot(z_mt, x_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    ax[1][0].plot(t_mt, vx_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    ax[1][1].plot(t_mt, vy_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    ax[1][2].plot(t_mt, vz_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    ax[2][0].plot(t_mt, vt_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    ax[2][1].plot(t_mt, dist_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    ax[2][2].plot(0,0, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    
                    # ax1[0].plot(x_mt, y_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    # ax1[1].plot(t_mt, dist_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
                    # ax1[2].plot(0,0, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')

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
            ax[0][0].plot(x, y, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            ax[0][1].plot(y, z, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            ax[0][2].plot(z, x, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            ax[1][0].plot(t, vx, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            ax[1][1].plot(t, vy, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            ax[1][2].plot(t, vz, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            ax[2][0].plot(t,vt, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            ax[2][1].plot(t,dist, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            ax[2][2].plot(0,0, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            
            # ax1[0].plot(x, y, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            # ax1[1].plot(t,dist, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
            # ax1[2].plot(0,0, ls=ls, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)

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
                    
                    ax[0][0].plot(x_d, y_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    ax[0][1].plot(y_d, z_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    ax[0][2].plot(z_d, x_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    ax[1][0].plot(t, vx_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    ax[1][1].plot(t, vy_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    ax[1][2].plot(t, vz_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    ax[2][0].plot(t, vt_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    ax[2][1].plot(t, dist_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    ax[2][2].plot(0, 0, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    
                    # ax1[0].plot(x_d, y_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    # ax1[1].plot(t, dist_d, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    # ax1[2].plot(0, 0, ls=ls_d, lw=3, label=r'$\mathrm{\sigma/m_{\chi} = %i\ cm^2/g}$'%sigma)
                    
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

                    ax[0][0].plot(x_tri, y_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
                    ax[0][1].plot(y_tri, z_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')			
                    ax[0][2].plot(z_tri, x_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
                    ax[1][0].plot(t, vx_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Drag}$')
                    ax[1][1].plot(t, vy_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
                    ax[1][2].plot(t, vz_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
                    ax[2][0].plot(t, vt_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
                    ax[2][1].plot(t, dist_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
                    ax[2][2].plot(0,0, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')

            ####### PLOTTING #######

            ####### ORBITS #######

            ax[0][0].plot(0,0,'k*',markersize=12) # , label=r'$\mathrm{Host\ Center}$')
            ax[0][0].plot(x[0],y[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # , label=r'$\mathrm{Orbit\ Start}$')
            ax[0][0].set_xlabel(r'$\mathrm{x\ (R_{200m})}$')
            ax[0][0].set_ylabel(r'$\mathrm{y\ (R_{200m})}$')
            # ax[0][0].legend(loc='lower left',fontsize=15)
            # ax[0][0].set_xlim(-rmax,rmax)
            # ax[0][0].set_ylim(-rmax,rmax)

            ax1[0].plot(0,0,'k*',markersize=12) # , label=r'$\mathrm{Host\ Center}$')
            ax1[0].plot(x[0],y[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # , label=r'$\mathrm{Orbit\ Start}$')
            ax1[0].set_xlabel(r'$x\ (R_{\rm 200m})$')
            ax1[0].set_ylabel(r'$y\ (R_{\rm 200m})$')
            

            ax[0][1].plot(0, 0, 'k*', markersize=12) # , label=r'$\mathrm{Host\ Center}$')
            ax[0][1].plot(y[0], z[0], '*', color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
            ax[0][1].set_xlabel(r'$\mathrm{y\ (R_{200m})}$')
            ax[0][1].set_ylabel(r'$\mathrm{z\ (R_{200m})}$')
            # ax[0][1].legend(loc='lower left',fontsize=15)
            # ax[0][1].set_xlim(-rmax,rmax)
            # ax[0][1].set_ylim(-rmax,rmax)

            ax[0][2].plot(0,0,'k*', markersize=12) # , label=r'$\mathrm{Host\ Center}$')
            ax[0][2].plot(z[0],x[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # , label=r'$\mathrm{Orbit\ Start}$')
            ax[0][2].set_xlabel(r'$\mathrm{z\ (R_{200m})}$')
            ax[0][2].set_ylabel(r'$\mathrm{x\ (R_{200m})}$')
            # ax[0][2].legend(loc='lower left',fontsize=15)
            # ax[0][2].set_xlim(-rmax,rmax)
            # ax[0][2].set_ylim(-rmax,rmax)

            ####### VELOCITIES #######

            ax[1][0].plot(t[0],vx[0],'*',color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
            ax[1][0].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
            ax[1][0].set_ylabel(r'$\mathrm{vx\ (km/s)}$')
            # ax[1][0].legend(loc='lower left',fontsize=15)
            # ax[1][0].set_xlim(t.min()-1,t.max()+1)

            ax[1][1].plot(t[0],vy[0],'*',color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
            ax[1][1].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
            ax[1][1].set_ylabel(r'$\mathrm{vy\ (km/s)}$')
            # ax[1][1].legend(loc='lower left',fontsize=15)
            # ax[1][1].set_xlim(t.min()-1,t.max()+1)

            ax[1][2].plot(t[0],vz[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
            ax[1][2].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
            ax[1][2].set_ylabel(r'$\mathrm{vz\ (km/s)}$')
            # ax[1][2].legend(loc='lower left',fontsize=15)
            # ax[1][2].set_xlim(t.min()-1,t.max()+1)

            ####### RADIUS #######

            ax[2][0].plot(t[0],vt[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
            ax[2][0].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
            ax[2][0].set_ylabel(r'$\mathrm{v\ (km/s)}$')
            # ax[2][0].set_xlim(t.min()-1,t.max()+1)
            # ax[2][0].set_ylim(0,rmax)

            ####### VELOCITY MAGNITUDE #######

            ax[2][1].plot(t[0],dist[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
            ax[2][1].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
            ax[2][1].set_ylabel(r'$\mathrm{r\ (R_{200m})}$')
            # ax[2][1].set_xlim(t.min()-1,t.max()+1)
            # ax[2][1].set_ylim(0,rmax)
            
            ax1[1].plot(t[0],dist[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
            ax1[1].set_xlabel(r'$t\ {\rm ( Gyr)}$')
            ax1[1].set_ylabel(r'$r\ (R_{\rm 200m})$')

            ax[2][2].legend(loc='center',fontsize=30)
            ax[2][2].get_xaxis().set_visible(False)
            ax[2][2].get_yaxis().set_visible(False)
            ax[2][2].spines['bottom'].set_color('white')
            ax[2][2].spines['top'].set_color('white') 
            ax[2][2].spines['right'].set_color('white')
            ax[2][2].spines['left'].set_color('white')
            
            '''
            ax1[2].legend(loc='center',fontsize=30)
            ax1[2].get_xaxis().set_visible(False)
            ax1[2].get_yaxis().set_visible(False)
            ax1[2].spines['bottom'].set_color('white')
            ax1[2].spines['top'].set_color('white') 
            ax1[2].spines['right'].set_color('white')
            ax1[2].spines['left'].set_color('white')
            '''
            
            # fig.delaxes(ax[2][2])
            # plt.draw()
            # handles, labels = ax[0][0].get_legend_handles_labels()
            # plt.figlegend(handles,labels,(0.725,0.1))

            ax[0][1].set_title(r'$\mathrm{Host\ %i\ (s=%.2f,\ q=%.2f),\ Subhalo\ %i}$' %(host_idx, host.s, host.q, sub_idx),fontsize=30,y=1.08)

            plt.savefig(HOMEDIR + '/sidm_orbit_calculation/src/plots/%i_%.0e_%i.png'%(host_idx,dt,sub_idx))

plt.show()
