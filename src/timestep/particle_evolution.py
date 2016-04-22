import numpy as np

'''Governs a test particle position and momentum'''

pass_through_center = 0.

class ParticlePosition :
    '''Use this to compute new particle position'''
    def __init__(self,initial_position=None,dt=1.e14) :
        # what should be here? why is this a class? (!!)
        # record steps here? (!!)
        self.current_position = initial_position
        self.position_array = [(initial_position[0],initial_position[1])]
        print 'initial position = ', (initial_position[0]*np.cos(initial_position[1]),initial_position[0]*np.sin(initial_position[1]))
        self.dt = dt

    def update_step(self,momentum=None) :
        self.current_position[1]+=momentum[1]*self.dt/self.current_position[0]
        self.current_position[0]+=momentum[0]*self.dt
        if self.current_position[0] < 0:
            self.current_position[0] = np.abs(self.current_position[0])
            self.current_position[1]+=np.pi
            pass_through_center = 1.
        else:
            pass_through_center = 0.
        # self.current_position[1] = self.current_position[1]%(2.*np.pi)
        print 'position =  %10.5g %10.5g' % (self.current_position[0],self.current_position[1])
        print 'position =  %10.5g %10.5g' % (self.current_position[0],self.current_position[1])
        self.position_array.append((self.current_position[0],self.current_position[1]))
        return None

    def _compute_next_position(self) :
        # what is this? (!!)
        pass

class ParticleMomentum :
    '''Use this to compute new particle momentum'''
    def __init__(self,initial_momentum=None,dt=1.e14) :
        self.current_momentum = initial_momentum
        self.momentum_array = [initial_momentum]
        self.dt = dt

    def update_step(self,position=None,force=None) :
        if pass_through_center:
            self.current_momentum[0] = -(self.current_momentum[0])+force[0]*dt
        else:
            self.current_momentum[0]+=force[0]*self.dt            
        self.current_momentum[1]+=force[1]*self.dt/position[0]
        self.momentum_array.append((self.current_momentum[0],self.current_momentum[1]))
        return None

def polar_to_cartesian(position=None,polar=None):
    pass 

def cartesian_to_polar(position=None,cartesian=None):
    pass

