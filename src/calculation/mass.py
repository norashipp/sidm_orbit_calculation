from __future__ import division
from scipy.integrate import tplquad
from scipy.integrate import quad
import numpy as np
import time
from colossus.halo.mass_so import R_to_M

from sidm_orbit_calculation.src.potentials.density import *
from sidm_orbit_calculation.src.utils.constants import *
# from sidm_orbit_calculation.src.utils.geometry import *

'''
def calculate_triaxial_mass(host,a,b):
    t0 = time.time()
    def integrand(z, y, x):
        r2 = x*x + y*y + z*z
        if r2>b*b:
            return 0
        else:
            return triaxial_NFW_density(x, y, z, host)
            
    gfun = lambda x: a
    hfun = lambda x: b
    qfun = lambda x, y: a
    rfun = lambda x, y: b

    res = tplquad(integrand,a,b,gfun,hfun,qfun,rfun,epsabs=1,epsrel=1)[0]

    t1 = time.time()
    print 'time passed = %.2f seconds' %(t1-t0)
    return res

def triaxial_integrand(phi, theta, r, host):
    sinth = np.sin(theta)
    costh = np.cos(theta)
    sinph = np.sin(phi)
    cosph = np.cos(phi)
    r2 = r*r
    return triaxial_NFW_density(r*cosph*sinth, r*sinph*sinth, r*costh, host) * sinth*r2
'''
def calculate_triaxial_mass(host,a,b):
    t0 = time.time()
    # func = lambda phi, theta, r: triaxial_NFW_density(r*np.cos(phi)*np.sin(theta), r*np.sin(phi)*np.sin(theta), r*np.cos(theta), host) * np.sin(theta)*r*r
    def integrand(phi, theta, r):
        sinth = np.sin(theta)
        costh = np.cos(theta)
        sinph = np.sin(phi)
        cosph = np.cos(phi)
        r2 = r*r
        return triaxial_NFW_density(r*cosph*sinth, r*sinph*sinth, r*costh, host) * sinth*r2

    gfun = lambda r: 0
    hfun = lambda r: np.pi
    qfun = lambda r, theta: 0
    rfun = lambda r, theta: 2*np.pi

    epsabs = 1.
    epsrel = 1.

    # res = tplquad(func,a,b,gfun,hfun,qfun,rfun)[0]
    res = tplquad(integrand, a, b, gfun, hfun, qfun, rfun, epsabs=epsabs, epsrel=epsrel)[0]
    t1 = time.time()
    print 'time passed = %.2f seconds' %(t1-t0)
    return res

def calculate_spherical_mass(host,a,b):
    # R = r/host.R_s
    func = lambda r: host.rho_s / (r/host.R_s * (1+r/host.R_s)**2) * r**2
    res = quad(func,a,b)[0]
    mass = res*4*np.pi
    return mass

def calculate_point_mass(host,a=0,b=0): # not ideal, clean up later
    return host.M

'''
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
'''

mass_dict = {'point_mass':calculate_point_mass,'spherical_NFW':calculate_spherical_mass,'triaxial_NFW':calculate_triaxial_mass, 'triaxial_NFW_BT':calculate_triaxial_mass}