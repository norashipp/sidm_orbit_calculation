import cPickle
import sys
import glob 

from sidm_orbit_calculation.src.utils.setup import *

host_idx = int(sys.argv[1])

# sub_idx = int(sys.argv[2])
integrator = 'leapfrog'
potential = 'spherical_NFW'
dt = 4e-3

infiles = glob.glob(HOMEDIR+'sidm_orbit_calculation/src/output/%i_%s_%s_%.0e_*.dat' %(host_idx,integrator,potential,dt))
infiles.sort()

sub_idx = 0
max_len = 0
max_file = None
for infile in infiles:
	f = open(infile,'rb')
	data = cPickle.load(f)
	f.close()

	# times,positions,momenta,gravity,drag,density,energy,host_idx,host_radius = data
	# times,positions,momenta,gravity,drag,density,energy,host_idx,potential,host_radius = data
	# times, positions = data
	curr_len = len(data[0])
	if curr_len > max_len and sub_idx != 101:
		print 'subhalo %i' %sub_idx
		max_len = curr_len
		max_file = infile
		print max_len
		print data[0][1]
		print max_file

	sub_idx+=1

print max_file
print max_len