import numpy as np
from gridding import *

# SkyImage is a holder that has both RA and DEC increasing with array index
# This convention is necessary for the FFT step
# However, to display this image in the traditional sky convention (North up,
# East to the left), you must set the first array element to the lower left
# corner *and* flip the array along the RA axis: `fliplr(data)`

class SkyImage:
    def __init__(self, data, ra, dec, lams):
        if len(data.shape) < 3:
            self.data = np.reshape(data, (data.shape[0], data.shape[1], 1))
        else:
            self.data = data
        self.ra = ra
        self.dec = dec
        self.lams = np.array(lams)

# data visibility is stored as VV[index, l]
# ie there is no regular grid, just a list
class Visibility:
    def __init__(self, VV, uu, vv, wgts, lams):
        self.VV = VV
        self.uu = uu
        self.vv = vv
        self.wgts = wgts
        self.lams = np.array(lams)

# model visibilities aren't weighted
# they also are stored VV[j, i, l]
class Model_Visibility:
    def __init__(self, VV, uu, vv, lams):
        self.VV = VV
        self.uu = uu
        self.vv = vv
        self.lams = np.array(lams)

class GCF_holder:
    def __init__(self, index_arr, gcf_arr, w_arr, uu, vv):
        self.index_arr = index_arr
        self.gcf_arr = gcf_arr
        self.w_arr = w_arr
        self.uu = uu
        self.vv = vv
