import numpy as np

class HostHalo:

    def __init__(self,M,potential):
        self.M = M
        self.potential = potential
        self.c_vir = self.concentration()
        self.R_200 = self.virial_radius()
        self.R_s = self.scale_radius()
        self.v_200 = self.virial_velocity()
        self.rho_s = self.scale_density()

    def concentration(self):
        c_vir = 4.5 # typical fro M~1e14, from Benedikt's paper, fig 5
        return c_vir

    def virial_radius(self):
        radius = 1.e22 # meters, ~10 Mpc
        return radius

    def scale_radius(self):
        return self.R_200/self.c_vir

    def virial_velocity(self):
        velocity = 1.e6 # m/s
        return velocity

    def scale_density(self):
        return self.M/(4*np.pi*self.R_s**3*(np.log(1+self.c_vir)-self.c_vir/(1+self.c_vir)))
        