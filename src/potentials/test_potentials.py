import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import gala.integrate as gi
import gala.dynamics as gd
import gala.potential as gp
from gala.units import UnitSystem

from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.potentials.test_spherical_potentials import *
from sidm_orbit_calculation.src.utils.constants import *

usys = UnitSystem(u.meter, u.second, u.kilogram, u.degree, u.meter/u.second)

host = HostHalo(M=1e14,potential='triaxial_NFW')
pot = gp.LeeSutoTriaxialNFWPotential(v_c = host.v*u.m/u.s, r_s = host.R_s*u.m, a = host.R*u.m, b = host.R*host.q*u.m, c = host.R*host.s*u.m, units=usys)

ics = gd.CartesianPhaseSpacePosition(pos=[1e22,0,0.]*u.m, vel=[0,3e5,0]*u.m/u.s)

orbit = pot.integrate_orbit(ics, dt=4e6*u.year, n_steps=5000, Integrator=gi.LeapfrogIntegrator)

print orbit.pos.shape

plt.figure()
plt.plot(orbit.pos[0],orbit.pos[1],'c-',linewidth=2,label='gala triaxial')
plt.xlabel('x')
plt.ylabel('y')
plt.show()