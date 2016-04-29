from __future__ import division
import numpy as np

class GetDragForce :
    def __init__(self, host=None) :
        self.host = host

        self.force = 0.
        self.force_array = []

    def calculate_drag_force(self, position=None, momentum=None) :
        sigma_mDM = 1.e-4 # cm2/g
        rho = self.host.density(position)
        return 0.25 * sigma_mDM * momentum**2 * rho