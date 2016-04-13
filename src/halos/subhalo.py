import numpy as np

class Subhalo:

    def __init__(self,mass_ratio=None,initial_position=None,initial_momentum=None):
        self.mass_ratio = mass_ratio

        self.position = initial_position
        self.momentum = initial_momentum