import matplotlib.pyplot as plt

from sidm_orbit_calculation.src.orbit_parameters.orbit_distributions import *
from sidm_orbit_calculation.src.halos.subhalo import *
from sidm_orbit_calculation.src.halos.host_halo import *
from sidm_orbit_calculation.src.utils.constants import *

host = HostHalo(1e14,'spherical_NFW')
subhalo = Subhalo(host,1e12,np.zeros(3),np.zeros(3))

total_inverse_cdf, radial_inverse_cdf = get_inverse_cdf(subhalo.host.M/M_sol, subhalo.M/M_sol)

vtot = []
vr = []
vth = []

i = 0
n = 10000
while i < n:
    x_total, x_radial = uniform(0,1,2)

    total_ratio = float(total_inverse_cdf(x_total))
    radial_ratio = float(radial_inverse_cdf(x_radial))

    theta_ratio = float(total_ratio * np.sqrt(1 - radial_ratio ** 2))

    vtot.append(total_ratio)
    vr.append(radial_ratio*total_ratio)
    vth.append(theta_ratio)

    i+=1
   
plt.figure()
plt.hist(vr,normed=True,bins=20)
plt.xlabel('v_r/v_200')

plt.figure()
plt.hist(vth,normed=True,bins=20)
plt.xlabel('v_theta/v200')

plt.figure()
plt.hist(vtot,normed=True,bins=20)
plt.xlabel('v_tot/v200')

plt.show()
