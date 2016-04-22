from __future__ import division
import pylab
import numpy as np

from sidm_orbit_calculation.src.utils.constants import *

class Plotting:

	def __init__(self,times=None,positions=None,momenta=None,forces=None,host=None):
		self.host = host
		self.times = times
		self.r = positions[:,0]
		self.phi = positions[:,1]%(2*np.pi)
		self.momenta = momenta
		self.forces = forces

	def orbit_color(self):
		x = self.r*np.cos(self.phi)*m_to_kpc
		y = self.r*np.sin(self.phi)*m_to_kpc
		__dt = self.times[1]-self.times[0]
		color = [tt/self.times.max() for tt in self.times]
		clr = pylab.cm.jet(color)
		pylab.figure()
		pylab.scatter(x,y,c=clr)
		pylab.plot(x[0],y[0],'k*',markersize=10)
		pylab.xlabel('x (kpc)')
		pylab.ylabel('y (kpc)')
		# pylab.colorbar()
		pylab.show()

	def orbit(self):
		x = self.r*np.cos(self.phi)
		y = self.r*np.sin(self.phi)
		pylab.figure()
		pylab.plot(x,y,'g-',markersize=10,linewidth=2)
		pylab.plot(x[0],y[0],'k*',markersize=10)
		pylab.xlabel('x')
		pylab.ylabel('y')
		# pylab.xlim([-1e22,1e22])
		pylab.show()

	def radial_position(self):
		pylab.figure()
		pylab.plot(self.times,self.r/self.host.R_200,'r-',markersize=10,linewidth=2)
		pylab.xlabel('time')
		pylab.ylabel('raidus/R200')
		pylab.show()

	def radial_position_color(self):
		# not working
		pylab.figure()
		pylab.scatter(self.times,self.r/self.host.R_200,c = pylab.cm.jet(np.log(self.times)/np.log(max(self.times))))
		pylab.xlabel('time')
		pylab.ylabel('raidus/R200')
		pylab.show()

	def angular_position(self):
		pylab.figure()
		pylab.plot(self.times,self.phi,'b',markersize=10,linewidth=2)
		pylab.xlabel('time')
		pylab.ylabel('phi')
		pylab.show()

	def gravitational_force(self):
		pylab.figure()
		pylab.plot(self.r/self.host.R_200,self.forces,'k.',markersize=10,linewidth=2)
		pylab.xlabel('radius/R200')
		pylab.ylabel('gravitational force')
		pylab.show()

	def radial_velocity(self):
		pylab.figure()
		pylab.plot(self.times,self.momenta[:,0]/self.host.v_200,'r',markersize=10,linewidth=2)
		pylab.xlabel('time')
		pylab.ylabel('velocity/v200')
		pylab.show()

