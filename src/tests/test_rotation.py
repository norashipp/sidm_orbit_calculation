import sys

from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.utils.geometry import *

host_idx = 40
potential = 'spherical_NFW'
host = HostHalo(host_idx,potential,subs=False)
host.update(host.cosmo.age(0))
gravity = GetGravitationalForce(host)

axis = np.array([sys.argv[1],sys.argv[2],sys.argv[3]],dtype=float)
norm = linalg.norm(axis)
host.ax = axis[0]/norm
host.ay = axis[1]/norm
host.az = axis[2]/norm

pr = rotate([-1,-1,-1],axis)

print '%.2f %.2f %.2f' %(pr[0], pr[1], pr[2])
