import numpy as np
import cPickle
import sys
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology

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
print infile
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

sub = subs[sub_idx]

mt_t = cosmo.age(1/sub.a-1)
# print mt_t[0],mt_t[-1]

mt_x = sub.rel_x
mt_y = sub.rel_y
mt_z = sub.rel_z

mt_vx = sub.rel_vx*1000*m_to_Mpc/s_to_Gyr
mt_vy = sub.rel_vy*1000*m_to_Mpc/s_to_Gyr
mt_vz = sub.rel_vz*1000*m_to_Mpc/s_to_Gyr

mt_dist = np.sqrt(mt_x**2 + mt_y**2 + mt_z**2)

# print mt_x[0], mt_y[0]
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

	dist_drag = np.sqrt(x_d**2 + y_d**2 + z_d**2)

	vx_d = velocities_drag[:,0]
	vy_d = velocities_drag[:,1]
	vz_d = velocities_drag[:,2]

triaxial = 1
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

gala = 1
if gala:
	f = open(HOMEDIR + 'sidm_orbit_calculation/src/output/gala_orbit_%i_%i.dat'%(host_idx,sub_idx))
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

fig, ax = plt.subplots(3,3,figsize=(20,20))
# plt.title(r'$\mathrm{Host\ %i,\ Subhalo\ %i}$' %(host_idx, sub_idx))

####### ORBITS #######

ax[0][0].plot(x, y, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
ax[0][0].plot(mt_x, mt_y, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[0][0].plot(0,0,'^k',markersize=12,label=r'$\mathrm{Host\ Center}$')
ax[0][0].plot(x[0],y[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if drag: ax[0][0].plot(x_d, y_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
if triaxial: ax[0][0].plot(x_tri, y_tri, 'k--', lw=3, label=r'$\mathrm{Triaxial}$')
if gala: ax[0][0].plot(x_g, y_g, '--', color='violet', lw=3, label=r'$\mathrm{Gala}$')
ax[0][0].set_xlabel(r'$\mathrm{x\ (Mpc)}$')
ax[0][0].set_ylabel(r'$\mathrm{y\ (Mpc)}$')
# ax[0][0].legend(loc='lower left',fontsize=15)
ax[0][0].set_xlim([-2,2])
ax[0][0].set_ylim([-2,2])

ax[0][1].plot(y, z, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
ax[0][1].plot(mt_y, mt_z, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[0][1].plot(0,0,'^k',markersize=12,label=r'$\mathrm{Host\ Center}$')
ax[0][1].plot(y[0],z[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if drag: ax[0][1].plot(y_d, z_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
if triaxial: ax[0][1].plot(y_tri, z_tri, 'k--', lw=3, label=r'$\mathrm{Triaxial}$')
if gala: ax[0][1].plot(y_g, z_g, '--', color='violet', lw=3, label=r'$\mathrm{Gala}$')
ax[0][1].set_xlabel(r'$\mathrm{y\ (Mpc)}$')
ax[0][1].set_ylabel(r'$\mathrm{z\ (Mpc)}$')
# ax[0][1].legend(loc='lower left',fontsize=15)
ax[0][1].set_xlim([-2,2])
ax[0][1].set_ylim([-2,2])

ax[0][2].plot(z, x, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
ax[0][2].plot(mt_z, mt_x, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[0][2].plot(0,0,'^k',markersize=12,label=r'$\mathrm{Host\ Center}$')
ax[0][2].plot(z[0],x[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if drag: ax[0][2].plot(z_d, x_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
if triaxial: ax[0][2].plot(z_tri, x_tri, 'k--', lw=3, label=r'$\mathrm{Triaxial}$')
if gala: ax[0][2].plot(z_g, x_g, '--', color='violet', lw=3, label=r'$\mathrm{Gala}$')
ax[0][2].set_xlabel(r'$\mathrm{z\ (Mpc)}$')
ax[0][2].set_ylabel(r'$\mathrm{x\ (Mpc)}$')
# ax[0][2].legend(loc='lower left',fontsize=15)
ax[0][2].set_xlim([-2,2])
ax[0][2].set_ylim([-2,2])

####### VELOCITIES #######

ax[1][0].plot(t, vx, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
ax[1][0].plot(mt_t, mt_vx, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[1][0].plot(t[0],vx[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if drag: ax[1][0].plot(t, vx_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
if triaxial: ax[1][0].plot(t, vx_tri, 'k--', lw=3, label=r'$\mathrm{Drag}$')
if gala: ax[1][0].plot(t_g, vx_g, '--', color='violet', lw=3, label=r'$\mathrm{Gala}$')
ax[1][0].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
ax[1][0].set_ylabel(r'$\mathrm{vx\ (Mpc/Gyr)}$')
# ax[1][0].legend(loc='lower left',fontsize=15)

ax[1][1].plot(t, vy, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
ax[1][1].plot(mt_t, mt_vy, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[1][1].plot(t[0],vy[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if drag: ax[1][1].plot(t, vy_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
if triaxial: ax[1][1].plot(t, vy_tri, 'k--', lw=3, label=r'$\mathrm{Triaxial}$')
if gala: ax[1][1].plot(t_g, vy_g, '--', color='violet', lw=3, label=r'$\mathrm{Gala}$')
ax[1][1].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
ax[1][1].set_ylabel(r'$\mathrm{vy\ (Mpc/Gyr)}$')
# ax[1][1].legend(loc='lower left',fontsize=15)

ax[1][2].plot(t, vz, 'c', lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
ax[1][2].plot(mt_t, mt_vz, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[1][2].plot(t[0],vz[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if drag: ax[1][2].plot(t, vz_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
if triaxial: ax[1][2].plot(t, vz_tri, 'k--', lw=3, label=r'$\mathrm{Triaxial}$')
if gala: ax[1][2].plot(t_g, vz_g, '--', color='violet', lw=3, label=r'$\mathrm{Gala}$')
ax[1][2].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
ax[1][2].set_ylabel(r'$\mathrm{vz\ (Mpc/Gyr)}$')
# ax[1][2].legend(loc='lower left',fontsize=15)

####### RADIUS #######

ax[2][1].plot(t,dist,'c',lw=3, label=r'$\mathrm{Orbit\ Calculation}$')
ax[2][1].plot(mt_t, mt_dist, 'g', lw=3, label=r'$\mathrm{Merger\ Tree}$')
ax[2][1].plot(t[0],dist[0],'c*',markersize=15,label=r'$\mathrm{Orbit\ Start}$')
if drag: ax[2][1].plot(t, dist_d, 'b--', lw=3, label=r'$\mathrm{Drag}$')
if triaxial: ax[2][1].plot(t, dist_tri, 'k--', lw=3, label=r'$\mathrm{Triaxial}$')
if gala: ax[2][1].plot(t_g, dist_g, '--', color='violet', lw=3, label=r'$\mathrm{Gala}$')
ax[2][1].set_xlabel(r'$\mathrm{t\ (Gyr)}$')
ax[2][1].set_ylabel(r'$\mathrm{r\ (Mpc)}$')

#plt.show()
plt.savefig(HOMEDIR + '/sidm_orbit_calculation/src/plots/%i_%s_%s_%.0e_%i.png'%(host_idx,potential,integrator,dt,sub_idx))
