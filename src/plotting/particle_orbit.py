import numpy as np
import sys
import matplotlib.pyplot as plt
# from colossus.cosmology import cosmology

from sidm_orbit_calculation.src.utils.setup import *
# from sidm_orbit_calculation.src.halos.host_halo import *

# my_cosmo = {'flat': True, 'H0': 70.0, 'Om0': 0.27,
#             'Ob0': 0.045714, 'sigma8': 0.82, 'ns': 0.96}
# cosmo = cosmology.setCosmology('my_cosmo', my_cosmo)

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
# particle_idx_array = np.array(sys.argv[2:],dtype=int)
# particle_idx_array = np.arange(0, 10)
particle_idx_array = [0]

potential = 'spherical_NFW'
# host = HostHalo(host_idx, potential)

particle_file = DATADIR + 'subs/particle_%i.dat' % host_idx
particles = np.loadtxt(particle_file, unpack=True)

for particle_idx in particle_idx_array:
    pi = particles[0]
    px = pi[::3]
    py = pi[1::3]
    pz = pi[2::3]

    fig, ax = plt.subplots(1, 3, figsize={10,3})
    ax[0].plot(px, py, lw=3, color='slateblue')
    ax[1].plot(py, pz, lw=3, color='slateblue')
    ax[2].plot(pz, px, lw=3, color='slateblue')

    ax[0].set_xlabel(r'$\mathrm{x}$')
    ax[1].set_xlabel(r'$\mathrm{y}$')
    ax[2].set_xlabel(r'$\mathrm{z}$')

    ax[0].set_ylabel(r'$\mathrm{y}$')
    ax[1].set_ylabel(r'$\mathrm{z}$')
    ax[2].set_ylabel(r'$\mathrm{x}$')

    plt.savefig('plots/particle_%i_%i.png' %(host_idx, particle_idx))
    # plt.show()
