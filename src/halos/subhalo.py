import numpy as np

from sidm_orbit_calculation.src.calculation.get_gravitational_force import *
from sidm_orbit_calculation.src.calculation.get_drag_force import *

class Subhalo:

    def __init__(self,host=None,mass_ratio=None,initial_position=None,initial_momentum=None):
        self.host = host

        self.mass_ratio = mass_ratio

        self.position = initial_position
        self.momentum = initial_momentum
        
        self.position_array = [(self.position[0],self.position[1])]
        self.momentum_array = [(self.momentum[0],self.momentum[1])]

        self.gravity = GetGravitationalForce(self.host)
        self.drag = GetDragForce(self.host)

       