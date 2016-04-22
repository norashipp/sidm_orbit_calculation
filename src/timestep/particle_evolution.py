import numpy as np

'''Governs a test particle position and momentum'''

'''EULER METHOD'''

# pass_through_center = 0.

class ParticlePosition :
    '''Use this to compute new particle position'''
    def __init__(self,initial_position=None,dt=1.e14):
        # what should be here? why is this a class? (!!)
        # record steps here? (!!)
        self.current_position = initial_position
        self.position_array = [(initial_position[0],initial_position[1])]
        print 'initial position = %10.5g %10.5g' % (initial_position[0]*np.cos(initial_position[1]),initial_position[0]*np.sin(initial_position[1]))
        self.dt = dt

    def update_step(self,momentum=None):
        self.current_position = self._compute_next_position(polar_to_cartesian(polar=self.current_position,position=self.current_position),momentum)
        # print 'position =  %10.5g %10.5g' % (self.current_position[0],self.current_position[1])
        self.position_array.append((self.current_position[0],self.current_position[1]))
        return None

    def _compute_next_position(self,position,momentum):
        position[1]+=momentum[1]*self.dt
        position[0]+=momentum[0]*self.dt
        return cartesian_to_polar(position)

class ParticleMomentum :
    '''Use this to compute new particle momentum'''
    def __init__(self,initial_momentum=None,dt=1.e14):
        self.current_momentum = initial_momentum
        self.momentum_array = [initial_momentum]
        self.dt = dt

    def update_step(self,position=None,force=None):
        print 'cartesian = ', polar_to_cartesian(polar=force,position=position)
        self._compute_next_momentum(polar_to_cartesian(polar=force,position=position))
        self.momentum_array.append((self.current_momentum[0],self.current_momentum[1]))
        return None

    def _compute_next_momentum(self,force):
        self.current_momentum[1]+=force[1]*self.dt
        self.current_momentum[0]+=force[0]*self.dt
        # still cartesian
        return None

def polar_to_cartesian(polar=None,position=None):
    x = polar[0]*np.cos(position[1]) # only while force is completely radial (!!)
    y = polar[0]*np.sin(position[1])
    return [x, y]

def cartesian_to_polar(cartesian=None):
    r = np.sqrt(cartesian[0]**2+cartesian[1]**2)
    phi = np.arctan(cartesian[1]/cartesian[0])
    if cartesian[1] < 0: phi+=np.pi # this is only ok as long as this function is only used on position (!!)
    return [r, phi]

