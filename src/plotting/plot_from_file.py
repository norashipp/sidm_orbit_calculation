from __future__ import division
import pylab
import numpy as np
import cPickle

from sidm_orbit_calculation.src.plotting.make_plots import *

f = open('../data/pickle.dat','rb')
data = cPickle.load(f)
f.close()

(times,positions,momenta,forces,host) = data

plot = Plotting(times=times,positions=positions,momenta=momenta,forces=forces,host=host)
plot.orbit()
plot.orbit_color()
plot.radial_position()
plot.angular_position()
plot.gravitational_force()
plot.radial_velocity()