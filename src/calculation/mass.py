from scipy.integrate import tplquad
from scipy.integrate import quad
import numpy as np
from colossus.halo.mass_so import R_to_M

from sidm_orbit_calculation.src.potentials.density import *
from sidm_orbit_calculation.src.utils.constants import *

def calculate_triaxial_mass(host):
    func = lambda phi, theta, r: triaxial_NFW_density(r, theta, phi, host)*np.sin(theta)*r**2
    a = 0
    b = host.R_s
    print R_to_M(host.R*m_to_kpc,0,'vir')
    gfun = lambda r: 0
    hfun = lambda r: np.pi
    qfun = lambda r, theta: 0
    rfun = lambda r, theta: 2*np.pi
    res = tplquad(func,a,b,gfun,hfun,qfun,rfun)[0]
    return res

def calculate_spherical_mass(host):
    # R = r/host.R_s
    func = lambda r: host.rho_s/(r/host.R_s*(1+r/host.R_s)**2)*r**2
    res = quad(func,0,host.R_s)[0]
    mass = res*4*np.pi
    return mass

def triaxial_NFW_density(r,theta,phi,host):
    x = r * np.cos(phi) * np.sin(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(theta)
    x2 = x*x
    y2 = y*y
    z2 = z*z
    q2 = host.q*host.q
    s2 = host.s*host.s
    r = np.sqrt(x2+y2/q2+z2/s2)
    alpha = 1
    eta = 3
    return host.rho_s/(r/host.R_s**alpha*(1+r/host.R_s)**(eta-alpha))
