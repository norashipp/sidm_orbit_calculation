'''
Collection of spherically symmetric test potentials
phi = phi(r)
'''

def spherical_harmonic_oscilator(A=1, B=1, r_sample=None) :

    return A + B*r_sample**2

def point_mass_potential() :
    pass

def isochrone_potential() :
    pass
