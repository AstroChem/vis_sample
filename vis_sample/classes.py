import numpy as np
from gridding import *

# SkyImage is a holder that has both RA and DEC increasing with array index
# This convention is necessary for the FFT step
# However, to display this image in the traditional sky convention (North up,
# East to the left), you must set the first array element to the lower left
# corner *and* flip the array along the RA axis: `fliplr(data)`

class SkyImage:
    def __init__(self, data, ra, dec, freqs):
        if len(data.shape) == 2:
            self.data = np.reshape(data, (data.shape[0], data.shape[1], 1))
        else:
            self.data = data
        self.ra = ra
        self.dec = dec
        self.freqs = np.array(freqs)



# data visibility is stored as VV[index, f]
# ie there is no regular grid, just a list
# therefore uu.shape == vv.shape

class Visibility:
    def __init__(self, VV, uu, vv, wgts, freqs):
        self.VV = VV
        self.uu = uu
        self.vv = vv
        self.wgts = wgts
        self.freqs = np.array(freqs)



# model visibilities aren't weighted and are on a regular grid
# they also are stored VV[v, u, f]
# uu and vv can be different sizes

class ModelVisibility:
    def __init__(self, VV, uu, vv, freqs):
        self.VV = VV
        self.uu = uu
        self.vv = vv
        self.freqs = np.array(freqs)


# The GcfHolder allows us to cache the several values for future calls
# it holds the pixel indexes for the dataset, the gcf values, and the total weights
# additionally we hold the uu, vv values, so that we can avoid importing the data multiple times

class GcfHolder:
    def __init__(self, index_arr, gcf_arr, w_arr, uu, vv):
        self.index_arr = index_arr
        self.gcf_arr = gcf_arr
        self.w_arr = w_arr
        self.uu = uu
        self.vv = vv
