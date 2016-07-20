from __future__ import division
import pylab
import numpy as np
import cPickle
import sys

from sidm_orbit_calculation.src.plotting.make_plots import *
from sidm_orbit_calculation.src.utils.setup import *

homedir = home_directory()

fname = sys.argv[1]
f = open(homedir+'sidm_orbit_calculation/src/output/' + fname,'rb')
data = cPickle.load(f)
f.close()

times,positions,momenta,gravity,drag,density,host = data
print 'number of saved time steps = ', times.shape

plot = Plotting(times=times, positions=positions, momenta=momenta, gravity=gravity, drag=drag, density=density ,host=host)
plot.orbit()
# plot.orbit_color()
# plot.radial_position()
# plot.angular_position()
# plot.gravitational_force()
# plot.radial_velocity()
