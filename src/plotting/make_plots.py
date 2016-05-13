from __future__ import division
import pylab
import numpy as np

from sidm_orbit_calculation.src.utils.constants import *

class Plotting:

	def __init__(self,times=None,positions=None,momenta=None,gravity=None,drag=None,density=None,host=None):
		self.host = host
		self.times = times
		self.x = positions[:,0]*m_to_kpc
		self.y = positions[:,1]*m_to_kpc
		self.r = np.sqrt(self.x**2+self.y**2)
		self.phi = np.arctan(self.y/self.x)
		self.momenta = momenta
		self.gravity = gravity
		self.drag = drag
		self.rho = density/(M_sol*m_to_kpc**3)

	def orbit_color(self):
		dt = self.times[1]-self.times[0]
		color = [tt/self.times.max() for tt in self.times]
		clr = pylab.cm.jet(color)
		pylab.figure()
		pylab.scatter(self.x,self.y,c=clr)
		pylab.plot(self.x[0],self.y[0],'k*',markersize=10)
		pylab.xlabel('x (kpc)')
		pylab.ylabel('y (kpc)')
		# pylab.colorbar()
		pylab.show()

	def orbit(self):
		pylab.figure()
		pylab.plot(self.x,self.y,'g-',markersize=10,linewidth=2)
		pylab.plot(self.x[0],self.y[0],'k*',markersize=10)
		pylab.xlabel('x')
		pylab.ylabel('y')
		pylab.xlim([-1e22*m_to_kpc,1e22*m_to_kpc])
		pylab.ylim([-1e22*m_to_kpc,1e22*m_to_kpc])
		pylab.show()

	def radial_position(self):
		pylab.figure()
		pylab.plot(self.times,self.r,'r-',markersize=10,linewidth=2)
		pylab.xlabel('time')
		pylab.ylabel('raidus (kpc)')
		pylab.show()

	def radial_position_color(self):
		# not working
		pylab.figure()
		pylab.scatter(self.times,self.r,c = pylab.cm.jet(np.log(self.times)/np.log(max(self.times))))
		pylab.xlabel('time')
		pylab.ylabel('raidus (kpc)')
		pylab.show()

	def angular_position(self):
		pylab.figure()
		pylab.plot(self.times,self.phi,'b',markersize=10,linewidth=2)
		pylab.xlabel('time')
		pylab.ylabel('phi')
		pylab.show()

	def gravitational_force(self):
		pylab.figure()
		pylab.plot(self.r[:-1],self.gravity,'k',linewidth=2)
		pylab.xlabel('radius (kpc)')
		pylab.ylabel('gravitational force')
		pylab.show()

	def drag_force(self):
		pylab.figure()
		pylab.plot(self.r[:-1]/self.host.R_200,np.sqrt(self.drag[:,0]**2+self.drag[:,1]**2),'c',linewidth=2)
		pylab.xlabel('radius/r200')
		pylab.ylabel('drag force')
		pylab.show()

	def density(self):
		pylab.figure()
		pylab.plot(self.r,self.rho,'m',linewidth=2)
		# pylab.plot(1,0.1*1e9,'k^',markersize=10) # check line passes reasonable value
		pylab.xlabel('r200 (kpc)')
		pylab.ylabel('density (Msol/kpc)')
		pylab.xscale('log')
		pylab.yscale('log')
		pylab.show()

	def radial_velocity(self):
		pylab.figure()
		pylab.plot(self.times,self.momenta[:,0]/self.host.v_200,'r',markersize=10,linewidth=2)
		pylab.xlabel('time')
		pylab.ylabel('velocity/v200')
		pylab.show()