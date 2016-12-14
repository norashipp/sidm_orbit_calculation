import numpy as np
import collections
import time

L = 125

class AbstractReader(object):
    def __init__(self, fname):
        self.fields = np.loadtxt(fname, unpack=True)
        self.id = np.array(self.fields[0], dtype=int)
        self.snap = np.array(self.fields[1], dtype=int)
        
    def _calc_indices(self):
        self.idxs = [0]
        for i in xrange(len(self.snap) - 1):
            if self.snap[i] > self.snap[i + 1]: self.idxs.append(1 + i)
        self.idxs.append(len(self.snap))

    def __len__(self):
        return len(self.idxs) - 1
    
class SubHalos(AbstractReader):
    def __init__(self, fname):
        AbstractReader.__init__(self, fname)
        (self.a, self.v_max, self.m_200m,
         self.x, self.y, self.z,
         self.vx, self.vy, self.vz) = self.fields[2:]
        
        self._calc_indices()

    def _get_raw_item(self, i):
        if i < 0 or i >= len(self):
            raise ValueError("Halo index %d out of bounds." % (i - 1))
        lo, hi = self.idxs[i], self.idxs[i+1]

        h = collections.namedtuple(
            "SubHalo", ["id", "snap", "a", "v_max", "m_200m",
                        "x", "y", "z", "rel_x", "rel_y", "rel_z",
                        "vx", "vy", "vz", "rel_vx", "rel_vy", "rel_vz"],
        )

        h.id = self.id[lo:hi]
        h.snap = self.snap[lo:hi]
        h.a = self.a[lo:hi]
        h.v_max = self.v_max[lo:hi]
        h.m_200m = self.m_200m[lo:hi]
        h.x = self.x[lo:hi]
        h.y = self.y[lo:hi]
        h.z = self.z[lo:hi]
        h.vx = self.vx[lo:hi]
        h.vy = self.vy[lo:hi]
        h.vz = self.vz[lo:hi]

        return h

    def get_host(self):
        return self._get_raw_item(0)

    def _wrap_coord(self, sub, host):
        diff = sub - host
        diff[diff > L/2] -= L
        diff[diff < -L/2] += L
        return diff
    
    def __getitem__(self, i):
        sub = self._get_raw_item(i + 1)
        host = self.get_host()
        n = min(len(host.id), len(sub.id))
        sub.x, sub.y, sub.z = sub.x[-n:], sub.y[-n:], sub.z[-n:]
        sub.vx, sub.vy, sub.vz = sub.vx[-n:], sub.vy[-n:], sub.vz[-n:]
        
        sub.rel_x = self._wrap_coord(sub.x, host.x[-n:])
        sub.rel_y = self._wrap_coord(sub.y, host.y[-n:])
        sub.rel_z = self._wrap_coord(sub.z, host.z[-n:])
        
        sub.rel_vx = sub.vx - host.vx[-n:]
        sub.rel_vy = sub.vy - host.vy[-n:]
        sub.rel_vz = sub.vz - host.vz[-n:]

        sub.id = sub.id[-n:]
        sub.snap = sub.snap[-n:]
        sub.a = sub.a[-n:]
        sub.m_200m = sub.m_200m[-n:]
        sub.v_max = sub.v_max[-n:]
        
        return sub
        
class HostHalos(AbstractReader):
    def __init__(self, fname):
        AbstractReader.__init__(self, fname)
        (self.a, self.m_200m, self.r_s,
         self.b_to_a, self.c_to_a, self.ax,
         self.ay, self.az) = self.fields[2:]
        
        self._calc_indices()

    def __getitem__(self, i):
        if i < 0 or i >= len(self):
            raise ValueError("Halo index %d out of bounds." % i)
        lo, hi = self.idxs[i], self.idxs[i+1]

        h = collections.namedtuple(
            "Halo", ["id", "snap", "a", "m_200m",
                     "r_s", "b_to_a", "c_to_a",
                     "ax", "ay", "az"],
        )

        h.id = self.id[lo:hi]
        h.snap = self.snap[lo:hi]
        h.a = self.a[lo:hi]
        h.m_200m = self.m_200m[lo:hi]
        h.r_s = self.r_s[lo:hi]
        h.b_to_a = self.b_to_a[lo:hi]
        h.c_to_a = self.c_to_a[lo:hi]
        h.ax = self.ax[lo:hi]
        h.ay = self.ay[lo:hi]
        h.az = self.az[lo:hi]

        norm = np.sqrt(h.ax*h.ax + h.ay*h.ay + h.az*h.az)
        h.ax /= norm
        h.ay /= norm
        h.az /= norm
        
        return h
    
if __name__ == "__main__":
    print time.time()
    hs = HostHalos("/scratch/midway/mansfield/data/merger_tree/clusters.dat")
    print time.time()
    print "%d hosts in clusters.dat" % len(hs)

    #host_idx = 40
    #snap = 60
    host_idx = 102
    snap = -1

    print "Host %d properties at snapshot %d:" % (host_idx, snap)
    print "    h.id = %d" % hs[host_idx].id[snap]
    print "    h.a = %.3g" % hs[host_idx].a[snap]
    print "    h.m_200m = %.7g" % hs[host_idx].m_200m[snap]
    print "    h.r_s = %.3g" % hs[host_idx].r_s[snap]
    print "    h.b_to_a = %.3g" % hs[host_idx].b_to_a[snap]
    print "    h.c_to_a = %.3g" % hs[host_idx].c_to_a[snap]
    print ("    major axis: (h.az, h.ay, h.az) = (%.3g, %.3g %.3g)" %
           (hs[host_idx].ax[snap],
            hs[host_idx].ay[snap],
            hs[host_idx].az[snap]))
    print

    print time.time()
    subs = SubHalos("subs/sub_%d.dat" % host_idx)
    print time.time()
    
    #sub_idx = 10
    sub_idx = 0

    print "%d subhalos for host %d" % (len(subs), host_idx)
    
    print ("Subhalo %d of halo %d at snapshot %d:" %
           (sub_idx, host_idx, snap))
    print "    h.a = %.3g" % subs[sub_idx].a[snap]
    print "    h.v_max = %.3g" % subs[sub_idx].v_max[snap]
    print "    h.m_200m = %.3g" % subs[sub_idx].m_200m[snap]
    print "    h.x = %.4g h.rel_x = %.4g" % (subs[sub_idx].x[snap],
                                             subs[sub_idx].rel_x[snap])
    print "    h.y = %.4g h.rel_y = %.4g" % (subs[sub_idx].y[snap],
                                             subs[sub_idx].rel_y[snap])
    print "    h.z = %.4g h.rel_z = %.4g" % (subs[sub_idx].z[snap],
                                             subs[sub_idx].rel_z[snap])
    print "    h.vx = %.4g h.rel_vx = %.4g" % (subs[sub_idx].vx[snap],
                                               subs[sub_idx].rel_vx[snap])
    
