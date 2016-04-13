import pylab

class Plotting:

	def __init__(self,times=None,positions=None,momenta=None,forces=None,host=None):
		self.host = host
		self.times = times
		self.positions = positions
		self.momenta = momenta
		self.forces = forces

	def orbit(self):
		pylab.figure()
		pylab.plot(self.positions[:,0]/self.host.R_200,self.positions[:,1],'g.',markersize=10)
		pylab.xlabel('phi')
		pylab.ylabel('raidus/R200')
		pylab.show()

	def radial_position(self):
		pylab.figure()
		pylab.plot(self.times,self.positions[:,0]/self.host.R_200,'r.',markersize=10)
		pylab.xlabel('time')
		pylab.ylabel('raidus/R200')
		pylab.show()

	def angular_position(self):
		pylab.figure()
		pylab.plot(self.times,self.positions[:,1],'b.',markersize=10)
		pylab.xlabel('time')
		pylab.ylabel('phi')
		pylab.show()

	def gravitational_force(self):
		pylab.figure()
		pylab.plot(self.positions[:-1,0]/self.host.R_200,self.forces,'k.',markersize=10)
		pylab.xlabel('radius/R200')
		pylab.ylabel('gravitational force')
		pylab.show()

	def radial_velocity(self):
		pylab.figure()
		pylab.plot(self.times,self.momenta[:,0]/self.host.v_200,'r.',markersize=10)
		pylab.xlabel('time')
		pylab.ylabel('velocity/v200')
		pylab.show()

