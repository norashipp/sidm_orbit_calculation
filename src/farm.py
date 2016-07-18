#!/usr/bin/env python

import os
import subprocess
import time
import numpy as np

from sidm_orbit_calculation.src.utils.constants import *


homedir = 'home/norashipp/sidm_orbit_calculation/'

username = 'norashipp'
jobname = 'sidm'

host_halo_mass = 1e13
dt = 1e4
tmax = 5e10 
n_particles = 10
subhalo_mass_array = np.array([5e11,1e12,5e12,1e13])
drag = 1
if drag:
    integrator = 'dissipative'
else:
    integrator = 'leapfrog'

potential = 'spherical_NFW'

for subhalo_mass in subhalo_mass_array:
    for index in range(0,n_particles):
        batch = 'sbatch --account=kicp --partition=westmere --job-name=%s --output=log.out --error=log.err --mem-per-cpu=24000 ' %(jobname)
        command = 'sim.py %.2e %.2e %.2e %.2e %s %s %i' %(host_halo_mass, subhalo_mass, dt, tmax, integrator, potential, index)
        command_queue = batch + command
        print command_queue
        outfile = homedir + 'src/output/%.1e_%.1e_%.1e_%.1e_%s_%s_%i.dat' %(host_halo_mass, subhalo_mass, dt*seconds_to_years, tmax*seconds_to_years, integrator, potential, index)

        if os.path.exists(outfile):
            print 'Outfile %s already exists, skipping...'%(outfile)
            continue

        os.system(command_queue) # Submit to queue
