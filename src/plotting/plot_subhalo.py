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
'''
fontsize = 15
plt.rc('savefig', dpi=dpi)
plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.rc('xtick.major', pad=5)
plt.rc('xtick.minor', pad=5)
plt.rc('ytick.major', pad=5)
plt.rc('ytick.minor', pad=5)
'''

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

t = times

x = positions[:,0]
y = positions[:,1]
z = positions[:,2]

vx = velocities[:,0]
vy = velocities[:,1]
vz = velocities[:,2]

dist = np.sqrt(x**2 + y**2 + z**2)
vt = np.sqrt(vx**2 + vy**2 + vz**2)

sub = subs[sub_idx]

t_mt = cosmo.age(1/sub.a-1)
# print t_mt[0],t_mt[-1]

h = 0.7
x_mt = sub.rel_x/(h*(1/sub.a))
y_mt = sub.rel_y/(h*(1/sub.a))
z_mt = sub.rel_z/(h*(1/sub.a))

vx_mt = sub.rel_vx*1000*m_to_Mpc/s_to_Gyr
vy_mt = sub.rel_vy*1000*m_to_Mpc/s_to_Gyr
vz_mt = sub.rel_vz*1000*m_to_Mpc/s_to_Gyr

dist_mt = np.sqrt(x_mt**2 + y_mt**2 + z_mt**2)
vt_mt = np.sqrt(vx_mt**2 + vy_mt**2 + vz_mt**2)

# print x_mt[0], y_mt[0]
# print x[0], y[0]

drag = 0
if drag:
	integrator = 'dissipative'
	infile_drag = HOMEDIR+'sidm_orbit_calculation/src/output/%i_%s_%s_%.0e_%i.dat' %(host_idx,integrator,potential,dt,sub_idx)
	f = open(infile_drag,'rb')
	data = cPickle.load(f)
	f.close()
	# times,positions,momenta,gravity,drag,density,energy,host_idx,host_radius = data
	# times,positions,momenta,gravity,drag,density,energy,host_idx,potential,host_radius = data
	_, positions_drag, velocities_drag = data
	x_d = positions_drag[:,0]
	y_d = positions_drag[:,1]
	z_d = positions_drag[:,2]

	dist_d = np.sqrt(x_d**2 + y_d**2 + z_d**2)

	vx_d = velocities_drag[:,0]
	vy_d = velocities_drag[:,1]
	vz_d = velocities_drag[:,2]

	vt_d = np.sqrt(vx_d**2 + vy_d**2 + vz_d**2)

triaxial = 0
if triaxial:
	integrator = 'leapfrog'
	potential = 'triaxial_NFW'
	infile_triaxial = HOMEDIR+'sidm_orbit_calculation/src/output/%i_%s_%s_%.0e_%i.dat' %(host_idx,integrator,potential,dt,sub_idx)
	f = open(infile_triaxial,'rb')
	data = cPickle.load(f)
	f.close()
	# times,positions,momenta,gravity,drag,density,energy,host_idx,host_radius = data
	# times,positions,momenta,gravity,drag,density,energy,host_idx,potential,host_radius = data
	_, positions_tri, velocities_tri = data
	x_tri = positions_tri[:,0]
	y_tri = positions_tri[:,1]
	z_tri = positions_tri[:,2]

	dist_tri = np.sqrt(x_tri**2 + y_tri**2 + z_tri**2)

	vx_tri = velocities_tri[:,0]
	vy_tri = velocities_tri[:,1]
	vz_tri = velocities_tri[:,2]

	vt_tri = np.sqrt(vx_tri**2 + vy_tri**2 + vz_tri**2)

gala = 1
if gala:
	f = open(HOMEDIR + 'sidm_orbit_calculation/src/output/gala_orbit_%i_spherical_NFW_%.0e_%i.dat'%(host_idx,dt,sub_idx))
	orbit = cPickle.load(f)
	f.close()
	x_g = orbit.w()[0]
	y_g = orbit.w()[1]
	z_g = orbit.w()[2]

	vx_g = orbit.w()[3]
	vy_g = orbit.w()[4]
	vz_g = orbit.w()[5]

	t_g = t[:-1]

	dist_g = np.sqrt(x_g**2 + y_g**2 + z_g**2)
	vt_g = np.sqrt(vx_g**2 + vy_g**2 + vz_g**2)

gala_tri = 0
if gala_tri:
	f = open(HOMEDIR + 'sidm_orbit_calculation/src/output/gala_orbit_%i_triaxial_NFW_%.0e_%i.dat'%(host_idx,dt,sub_idx))
	orbit = cPickle.load(f)
	f.close()
	x_gt = orbit.w()[0]
	y_gt = orbit.w()[1]
	z_gt = orbit.w()[2]

	vx_gt = orbit.w()[3]
	vy_gt = orbit.w()[4]
	vz_gt = orbit.w()[5]

	t_gt = t[:-1]

	dist_gt = np.sqrt(x_gt**2 + y_gt**2 + z_gt**2)
	vt_gt = np.sqrt(vx_gt**2 + vy_gt**2 + vz_gt**2)

constant = 0
if constant:
	integrator = 'leapfrog'
	potential = 'spherical_NFW'
	infile = HOMEDIR+'sidm_orbit_calculation/src/output/%i_%s_%s_%.0e_%i_constant_potential.dat' %(host_idx,integrator,potential,dt,sub_idx)
	f = open(infile,'rb')
	data = cPickle.load(f)
	f.close()
	_, positions, velocities = data
	x_const = positions[:,0]
	y_const = positions[:,1]
	z_const = positions[:,2]

	dist_const = np.sqrt(x_const**2 + y_const**2 + z_const**2)

	vx_const = velocities[:,0]
	vy_const = velocities[:,1]
	vz_const = velocities[:,2]

	vt_const = np.sqrt(vx_const**2 + vy_const**2 + vz_const**2)



'''
# COLORS
NCURVES = 6
values = range(NCURVES)
viridis = cm = plt.get_cmap('Blues') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=viridis)

c = scalarMap.to_rgba(values[0])
c_mt = scalarMap.to_rgba(values[1])
c_d = scalarMap.to_rgba(values[2])
c_tri = scalarMap.to_rgba(values[3])
c_g = scalarMap.to_rgba(values[4])
c_gt  = scalarMap.to_rgba(values[5])
'''

c_mt = 'g'
c = 'r'
c_tri = c
c_const = c
c_d = c
c_g = 'b'
c_gt = c_g

ls_mt = '-'
ls = '--'
ls_tri = '--'
ls_const = '-.'
ls_d = '-.'
ls_g = '--'
ls_gt = '--'

fig, ax = plt.subplots(3,3,figsize=(20,20))
# plt.title(r'$\mathrm{Host\ %i,\ Subhalo\ %i}$' %(host_idx, sub_idx))

####### ORBITS #######

rmax = 2

ax[0][0].plot(x_mt, y_mt, color=c_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[0][0].plot(x, y, color=c, ls=ls, lw=3, label=r'$\mathrm{Spherical\ NFW}$')
if triaxial: ax[0][0].plot(x_tri, y_tri, color=c_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
if constant: ax[0][0].plot(x_const, y_const, color=c_const, ls=ls_const, lw=3, label=r'$\mathrm{Constant\ Potential}$')
if drag: ax[0][0].plot(x_d, y_d, color=c_d, ls=ls_d, lw=3, label=r'$\mathrm{Drag}$')
if gala: ax[0][0].plot(x_g, y_g, color=c_g, ls=ls_g, lw=3, label=r'$\mathrm{Gala\ Spherical}$')
if gala_tri: ax[0][0].plot(x_gt, y_gt, color=c_gt, ls=ls_gt, lw=3, label=r'$\mathrm{Gala\ Triaxial}$')
ax[0][0].plot(0,0,'k*',markersize=12) # , label=r'$\mathrm{Host\ Center}$')
ax[0][0].plot(x[0],y[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # , label=r'$\mathrm{Orbit\ Start}$')
ax[0][0].set_xlabel(r'$\mathrm{x\ (Mpc)}$')
ax[0][0].set_ylabel(r'$\mathrm{y\ (Mpc)}$')
# ax[0][0].legend(loc='lower left',fontsize=15)
ax[0][0].set_xlim(-rmax,rmax)
ax[0][0].set_ylim(-rmax,rmax)

ax[0][1].plot(y_mt, z_mt, color=c_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[0][1].plot(y, z, color=c, ls=ls, lw=3, label=r'$\mathrm{Spherical\ NFW}$')
if triaxial: ax[0][1].plot(y_tri, z_tri, color=c_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
if constant: ax[0][1].plot(y_const, z_const, color=c_const, ls=ls_const, lw=3, label=r'$\mathrm{Constant\ Potential}$')
if drag: ax[0][1].plot(y_d, z_d, color=c_d, ls=ls_d, lw=3, label=r'$\mathrm{Drag}$')
if gala: ax[0][1].plot(y_g, z_g, color=c_g, ls=ls_g, lw=3, label=r'$\mathrm{Gala\ Spherical}$')
if gala_tri: ax[0][1].plot(y_gt, z_gt, color=c_gt, ls=ls_gt, lw=3, label=r'$\mathrm{Gala\ Triaxial}$')
ax[0][1].plot(0, 0, 'k*', markersize=12) # , label=r'$\mathrm{Host\ Center}$')
ax[0][1].plot(y[0], z[0], '*', color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
ax[0][1].set_xlabel(r'$\mathrm{y\ (Mpc)}$')
ax[0][1].set_ylabel(r'$\mathrm{z\ (Mpc)}$')
# ax[0][1].legend(loc='lower left',fontsize=15)
ax[0][1].set_xlim(-rmax,rmax)
ax[0][1].set_ylim(-rmax,rmax)

ax[0][2].plot(z_mt, x_mt, color=c_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[0][2].plot(z, x, color=c, ls=ls, lw=3, label=r'$\mathrm{Spherical\ NFW}$')
if triaxial: ax[0][2].plot(z_tri, x_tri, color=c_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
if constant: ax[0][2].plot(z_const, x_const, color=c_const, ls=ls_const, lw=3, label=r'$\mathrm{Constant\ Potential}$')
if drag: ax[0][2].plot(z_d, x_d,color=c_d, ls=ls_d, lw=3, label=r'$\mathrm{Drag}$')
if gala: ax[0][2].plot(z_g, x_g, color=c_g, ls=ls_g, lw=3, label=r'$\mathrm{Gala\ Spherical}$')
if gala_tri: ax[0][2].plot(z_gt, x_gt, color=c_gt, ls=ls_gt, lw=3, label=r'$\mathrm{Gala\ Triaxial}$')
ax[0][2].plot(0,0,'k*', markersize=12) # , label=r'$\mathrm{Host\ Center}$')
ax[0][2].plot(z[0],x[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # , label=r'$\mathrm{Orbit\ Start}$')
ax[0][2].set_xlabel(r'$\mathrm{z\ (Mpc)}$')
ax[0][2].set_ylabel(r'$\mathrm{x\ (Mpc)}$')
# ax[0][2].legend(loc='lower left',fontsize=15)
ax[0][2].set_xlim(-rmax,rmax)
ax[0][2].set_ylim(-rmax,rmax)

####### VELOCITIES #######

ax[1][0].plot(t_mt, vx_mt, color=c_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[1][0].plot(t, vx, color=c, ls=ls, lw=3, label=r'$\mathrm{Spherical\ NFW}$')
if triaxial: ax[1][0].plot(t, vx_tri, color=c_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Drag}$')
if constant: ax[1][0].plot(t, vx_const, color=c_const, ls=ls_const, lw=3, label=r'$\mathrm{Constant\ Potential}$')
if drag: ax[1][0].plot(t, vx_d,color=c_d, ls=ls_d, lw=3, label=r'$\mathrm{Drag}$')
if gala: ax[1][0].plot(t_g, vx_g, color=c_g, ls=ls_g, lw=3, label=r'$\mathrm{Gala\ Spherical}$')
if gala_tri: ax[1][0].plot(t_gt, vx_gt, color=c_gt, ls=ls_gt, lw=3, label=r'$\mathrm{Gala\ Triaxial}$')
ax[1][0].plot(t[0],vx[0],'*',color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
ax[1][0].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
ax[1][0].set_ylabel(r'$\mathrm{vx\ (Mpc/Gyr)}$')
# ax[1][0].legend(loc='lower left',fontsize=15)
ax[1][0].set_xlim(t.min()-1,t.max()+1)

ax[1][1].plot(t_mt, vy_mt, color=c_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[1][1].plot(t, vy, color=c, ls=ls, lw=3, label=r'$\mathrm{Spherical\ NFW}$')
if triaxial: ax[1][1].plot(t, vy_tri, color=c_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
if constant: ax[1][1].plot(t, vy_const, color=c_const, ls=ls_const, lw=3, label=r'$\mathrm{Constant\ Potential}$')
if drag: ax[1][1].plot(t, vy_d, color=c_d, ls=ls_d, lw=3, label=r'$\mathrm{Drag}$')
if gala: ax[1][1].plot(t_g, vy_g, color=c_g, ls=ls_g, lw=3, label=r'$\mathrm{Gala\ Spherical}$')
if gala_tri: ax[1][1].plot(t_gt, vy_gt, color=c_gt, ls=ls_gt, lw=3, label=r'$\mathrm{Gala\ Triaxial}$')
ax[1][1].plot(t[0],vy[0],'*',color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
ax[1][1].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
ax[1][1].set_ylabel(r'$\mathrm{vy\ (Mpc/Gyr)}$')
# ax[1][1].legend(loc='lower left',fontsize=15)
ax[1][1].set_xlim(t.min()-1,t.max()+1)

ax[1][2].plot(t_mt, vz_mt, color=c_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[1][2].plot(t, vz, color=c, ls=ls, lw=3, label=r'$\mathrm{Spherical\ NFW}$')
if triaxial: ax[1][2].plot(t, vz_tri, color=c_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
if constant: ax[1][2].plot(t, vz_const, color=c_const, ls=ls_const, lw=3, label=r'$\mathrm{Constant\ Potential}$')
if drag: ax[1][2].plot(t, vz_d, color=c_d, ls=ls_d, lw=3, label=r'$\mathrm{Drag}$')
if gala: ax[1][2].plot(t_g, vz_g, color=c_g, ls=ls_g, lw=3, label=r'$\mathrm{Gala\ Spherical}$')
if gala_tri: ax[1][2].plot(t_gt, vz_gt, color=c_gt, ls=ls_gt, lw=3, label=r'$\mathrm{Gala\ Triaxial}$')
ax[1][2].plot(t[0],vz[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
ax[1][2].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
ax[1][2].set_ylabel(r'$\mathrm{vz\ (Mpc/Gyr)}$')
# ax[1][2].legend(loc='lower left',fontsize=15)
ax[1][2].set_xlim(t.min()-1,t.max()+1)

####### RADIUS #######

ax[2][0].plot(t_mt, vt_mt, color=c_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[2][0].plot(t,vt, color=c, ls=ls, lw=3, label=r'$\mathrm{Spherical\ NFW}$')
if triaxial: ax[2][0].plot(t, vt_tri, color=c_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
if constant: ax[2][0].plot(t, vt_const, color=c_const, ls=ls_const, lw=3, label=r'$\mathrm{Constant\ Potential}$')
if drag: ax[2][0].plot(t, vt_d, color=c_d, ls=ls_d, lw=3, label=r'$\mathrm{Drag}$')
if gala: ax[2][0].plot(t_g, vt_g, color=c_g, ls=ls_g, lw=3, label=r'$\mathrm{Gala\ Spherical}$')
if gala_tri: ax[2][0].plot(t_gt, vt_gt, color=c_gt, ls=ls_gt, lw=3, label=r'$\mathrm{Gala\ Triaxial}$')
ax[2][0].plot(t[0],vt[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
ax[2][0].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
ax[2][0].set_ylabel(r'$\mathrm{v\ (Mpc/Gyr)}$')
ax[2][0].set_xlim(t.min()-1,t.max()+1)
# ax[2][0].set_ylim(0,rmax)

####### RADIUS #######

ax[2][1].plot(t_mt, dist_mt, color=c_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[2][1].plot(t,dist, color=c, ls=ls, lw=3, label=r'$\mathrm{Spherical\ NFW}$')
if triaxial: ax[2][1].plot(t, dist_tri, color=c_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
if constant: ax[2][1].plot(t, dist_const, color=c_const, ls=ls_const, lw=3, label=r'$\mathrm{Constant\ Potential}$')
if drag: ax[2][1].plot(t, dist_d, color=c_d, ls=ls_d, lw=3, label=r'$\mathrm{Drag}$')
if gala: ax[2][1].plot(t_g, dist_g, color=c_g, ls=ls_g, lw=3, label=r'$\mathrm{Gala\ Spherical}$')
if gala_tri: ax[2][1].plot(t_gt, dist_gt, color=c_gt, ls=ls_gt, lw=3, label=r'$\mathrm{Gala\ Triaxial}$')
ax[2][1].plot(t[0],dist[0],'*', color=c, ls=ls, markeredgecolor=None, markersize=15) # ,label=r'$\mathrm{Orbit\ Start}$')
ax[2][1].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
ax[2][1].set_ylabel(r'$\mathrm{r\ (Mpc)}$')
ax[2][1].set_xlim(t.min()-1,t.max()+1)
ax[2][1].set_ylim(0,rmax)

ax[2][2].plot(0,0, color=c_mt, ls=ls_mt, lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[2][2].plot(0,0, color=c, ls=ls, lw=3, label=r'$\mathrm{Spherical\ NFW}$')
if triaxial: ax[2][2].plot(0,0, color=c_tri, ls=ls_tri, lw=3, label=r'$\mathrm{Triaxial\ NFW}$')
if constant: ax[2][2].plot(0,0, color=c_const, ls=ls_const, lw=3, label=r'$\mathrm{Constant\ Potential}$')
if drag: ax[2][2].plot(0,0, color=c_d, ls=ls_d, lw=3, label=r'$\mathrm{Drag}$')
if gala: ax[2][2].plot(0,0, color=c_g, ls=ls_g, lw=3, label=r'$\mathrm{Gala\ Spherical}$')
if gala_tri: ax[2][2].plot(0,0, color=c_gt, ls=ls_gt, lw=3, label=r'$\mathrm{Gala\ Triaxial}$')

ax[2][2].legend(loc='lower center',fontsize=30)
ax[0][1].set_title(r'$\mathrm{Host\ %i,\ Subhalo\ %i}$' %(host_idx,sub_idx),fontsize=30)

plt.show()
# plt.savefig(HOMEDIR + '/sidm_orbit_calculation/src/plots/%i_%s_%.0e_%i.png'%(host_idx,integrator,dt,sub_idx))
