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

times,positions,momenta,gravity,drag,density,energy,host = data
print 'number of saved time steps = %i' % positions.shape[0]
print 'initial position = %.3g, %.3g, %.3g' % (positions[0][0],positions[0][1],positions[0][2])
print 'initial velocity = %.3g, %.3g, %.3g' % (momenta[0][0],momenta[0][1],momenta[0][2])

# print 'initial gravitational force = %.3g'  %gravity[0]

radius = np.sqrt(np.sum(positions**2,axis=1))
print radius.shape
print times.shape
print 'minimum radius = %.5e' %radius.min()
print 'maximum radius = %.5e' %radius.max()
print 'difference = %.5e' %(radius.max()-radius.min())

plot = Plotting(times=times, positions=positions, momenta=momenta, gravity=gravity, drag=drag, density=density, energy=energy, host=host)
plot.orbit()
plot.plot_energy()
# plot.y_position()
# plot.orbit_color()
plot.radius()
# plot.angular_position()
# plot.gravitational_force()
# plot.radial_velocity()
