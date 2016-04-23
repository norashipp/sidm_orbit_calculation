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
        return None

    def _compute_next_position(self,momentum):
        self.current_position[1]+=momentum[1]*self.dt
        self.current_position[0]+=momentum[0]*self.dt
        return None

class ParticleMomentum:
    '''Use this to compute new particle momentum'''
    def __init__(self,initial_momentum=None,gravity=None,dt=1.e14):
        self.current_momentum = initial_momentum
        self.momentum_array = [initial_momentum]
        self.dt = dt
        self.gravity = gravity # does this make sense?

    def update_step(self,position=None,force=None):
        force = self._update_force(position=position)
        self._compute_next_momentum(force)
        self.momentum_array.append((self.current_momentum[0],self.current_momentum[1]))
        return None

    def _compute_next_momentum(self,force):
        self.current_momentum[1]+=force[1]*self.dt
        self.current_momentum[0]+=force[0]*self.dt
        return None

    def _update_force(self,position=None):
        self.gravity.calculate_gravitational_force(position=position)
        phi = calculate_angle(position=position)
        force = self.gravity.vector(phi=phi)
        return force

def calculate_angle(position=None):
    return np.arctan2(position[1],position[0])

def polar_to_cartesian(polar=None,position=None):
    x = polar[0]*np.cos(position[1]) # only while force is completely radial (!!)
    y = polar[0]*np.sin(position[1])
    return [x, y]

def cartesian_to_polar(cartesian=None):
    r = np.sqrt(cartesian[0]**2+cartesian[1]**2)
    phi = np.arctan(cartesian[1]/cartesian[0])
    # if cartesian[1] < 0: phi+=np.pi # this is only ok as long as this function is only used on position (!!)
    return [r, phi]

