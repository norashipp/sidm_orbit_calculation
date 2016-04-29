import numpy as np

'''Governs a test particle position and momentum'''

'''EULER METHOD'''

# pass_through_center = 0.

class ParticlePosition:
    '''Use this to compute new particle position'''
    def __init__(self,initial_position=None,dt=1.e14):
        self.current_position = initial_position
        self.position_array = [(initial_position[0],initial_position[1])]
        # print 'initial position = %10.5g %10.5g' % (initial_position[0],initial_position[0])
        self.dt = dt

    def update_step(self,momentum=None):
        self._compute_next_position(momentum)
        self.position_array.append((self.current_position[0],self.current_position[1]))

    def _compute_next_position(self,momentum):
        self.current_position+=momentum*self.dt
        
class ParticleMomentum:
    '''Use this to compute new particle momentum'''
    def __init__(self,integration_method='leapfrog',initial_momentum=None,gravity=None,dt=1.e14,position=None):
        self.current_momentum = initial_momentum

        self.dt = dt
        self.gravity = gravity # does this make sense?
        
        if integration_method == 'leapfrog':
            self._first_momentum_step(position=position)
            print 'leapfrog'
    
        self.momentum_array = [(self.current_momentum[0],self.current_momentum[1])]
        
    def update_step(self,position=None):
        force = self._update_force(position=position)
        self._compute_next_momentum(force)
        self.momentum_array.append((self.current_momentum[0],self.current_momentum[1]))
        
    def _compute_next_momentum(self,force=None):
        self.current_momentum+=force*self.dt
        
    def _update_force(self,position=None):
        self.gravity.calculate_gravitational_force(position=position)
        phi = calculate_angle(position=position)
        force = self.gravity.vector(phi=phi)
        return force
        
    def _first_momentum_step(self,position=None):
        partial = self._partial_force(position=position)
        force = self._update_force(position=position)
        self.current_momentum+=self.dt/2*force+(self.dt/2)**2*partial
        
    def _partial_force(self,position=None):
        return self.gravity.calculate_partial_force(position=position)

def calculate_angle(position=None):
    return np.arctan2(position[1],position[0])
