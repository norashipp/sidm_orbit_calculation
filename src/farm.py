#!/usr/bin/env python

import os
import subprocess
import time
import numpy as np
import sys

from sidm_orbit_calculation.src.utils.constants import *
from sidm_orbit_calculation.src.utils.setup import *

username = 'norashipp'
jobname = 'sidm'

host_idx = int(sys.argv[1])
dt = 4e-3
drag = 0
if drag:
    integrator = 'dissipative'
else:
    integrator = 'leapfrog'
potential = 'triaxial_NFW'

batch = 'sbatch --account=kicp --partition=kicp --job-name=%s --output=log.out --error=log.err --mem-per-cpu=5000 ' %(jobname)
command = 'sim.py %i %.2e %s %s' %(host_idx, dt, integrator, potential)
command_queue = batch + command
print command_queue

'''
outfile = HOMEDIR + 'sidm_orbit_calculation/src/output/%i_%s_%s_%.2e_%i.dat' %(host_idx, integrator, potential, dt, sub_idx)

if os.path.exists(outfile):
    print 'Outfile %s already exists, skipping...'%(outfile)
    continue
'''

os.system(command_queue) # Submit to queue
