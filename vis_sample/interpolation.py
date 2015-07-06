# Here we do the interpolation of the fft'd image, and sample at each visibility point.
# Instead of analytically calculating the gcf for points around each (u,v) sample point
# we estimate the gcf by precalculating the gcf for a dense grid (1000x image res) and
# then use these gcf values (for a massive speed increase)

import numpy as np
from gridding import *
from classes import *
from numpy.lib.stride_tricks import as_strided



# Calculate the gcf values for 7 fft image points around the dense grid points
# This is calculated as a 1d problem, calculating u and v weights separately
# to be multiplied in the future as an outer product.
#
# dense_grid_gcf is calculated only once, and then gcf values are retrieved from it

def calc_dense_grid_gcf():
    dense_grid = np.linspace(-0.5, 0.5, 1000)
    dense_grid_gcf = np.zeros((1000,6))
    for i, grid_pnt in enumerate(dense_grid):
        eta = (np.arange(-3,3) + grid_pnt)/3.0
        dense_grid_gcf[i] = gcffun(eta)
    return dense_grid_gcf




# Caching the gcf values and data indices for a given dataset of (u,v) visibilities results
# in a ~20% speed increase.  This is implemented through the GCF_holder class.
#
# Here we calculated all the gcf values for the visibilities and then pass the values
# and indices to a GCF_holder

def create_gcf_holder(uu, vv, vis):
    nvis = vv.shape[0]
    viscnt = vis.uu.shape[0]
    dense_grid_gcf = calc_dense_grid_gcf()

    # 1. Find the nearest pixel points in the fft'd image, and push to the index array
    du = vis.uu[1] - vis.uu[0]
    dv = vis.vv[1] - vis.vv[0]
    iu0 = viscnt/2+np.round((uu/du)).astype(int)
    iv0 = viscnt/2+np.round((vv/dv)).astype(int)
    index_arr = np.transpose(np.array((iu0, iv0)))

    # 2. Find the relative distance to this point (should be -0.5 du/v < val < 0.5 du/v) 
    u0 = uu - vis.uu[iu0]
    v0 = vv - vis.vv[iv0]

    # 3. Now we take the relative distance and find the nearest dense grid point
    iu0_grid = 500+(1000*u0/du).astype(int)
    iv0_grid = 500+(1000*v0/dv).astype(int)

    # 4. Pull the gcf vals for the nearest 7 pixels around this dense grid point
    uw = dense_grid_gcf[iu0_grid,:]
    vw = dense_grid_gcf[iv0_grid,:]

    # 5. Calculate the weights (outer product) and an array of the total weights 
    all_weights = np.einsum('...a,...b->...ab', vw, uw)
    w_arr = np.einsum('...ab->...', all_weights)

    return GCF_holder(index_arr, all_weights, w_arr, uu, vv)




# The interpolation call will first calculate the gcf holder (if not provided), and then use
# those values to calculate the interpolated visibilities. 
#
# If the return_gcf flag is set, it returns the gcf holder (possibly useful for future interpolations 
# where the gcf calculation represents a significant fraction of the computation time).

def interpolate_uv(uu, vv, vis, gcf_holder='', return_gcf=False):
    # create gcf_holder if one isn't provided
    if (gcf_holder==""):
        t0 = time.time()
        gcf_holder = create_gcf_holder(uu, vv, vis)
        t1 = time.time()
        print "gcf creation time = " + str(t1-t0)

    # create the new interpolated visibility holder
    nvis = vv.shape[0]
    interp_vis = np.zeros((nvis, vis.lams.shape[0]), dtype='complex')
    
    # read out the pixel indices for windowing
    u0 = gcf_holder.index_arr[:, 0].astype(int)-3
    v0 = gcf_holder.index_arr[:, 1].astype(int)-3

    # iterate through the channels and multiply weights by the windowed pixels
    # note that we could theoretically vectorize by channel as well, but this
    # produces a memory error - maybe the array is too large?
    for l in range(vis.lams.shape[0]):
        VV_chan = vis.VV[:,:,l]
        VV = as_strided(VV_chan, shape=(VV_chan.shape[0]-5, VV_chan.shape[1]-5, 6, 6), strides=VV_chan.strides * 2)
        interp_vis[:,l] = np.einsum("...ab->...", gcf_holder.gcf_arr*VV[v0,u0])/gcf_holder.w_arr

    if (return_gcf==True):
        return interp_vis, gcf_holder
    else:
        return interp_vis
