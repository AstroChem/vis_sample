# Here we do the interpolation of the fft'd image, and sample at each visibility point.
# Instead of analytically calculating the gcf for points around each (u,v) sample point
# we estimate the gcf by precalculating the gcf for a dense grid (1000x image res) and
# then use these gcf values (for a massive speed increase)
import numpy as np
from gridding import *
from classes import *
from numpy.lib.stride_tricks import as_strided



# Calculate the gcf values for 5 fft image points around the dense grid points
# This is calculated as a 1d problem, calculating u and v weights separately
# to be multiplied in the future as an outer product.
#
# dense_grid_gcf is calculated only once, and then gcf values are retrieved from it

def calc_dense_grid_gcf():
    dense_grid = np.linspace(-0.5, 0.5, 1001)
    dense_grid_gcf = np.zeros((1001,5))
    for i, grid_pnt in enumerate(dense_grid):
        eta = (np.arange(-2,3) + grid_pnt)/2.5
        dense_grid_gcf[i] = gcffun(eta)
    return dense_grid_gcf




# Caching the gcf values and data indices for a given dataset of (u,v) visibilities results
# in a ~20% speed increase for each channel.  This is implemented through the GCF_holder class.
#
# Here we calculated all the gcf values for the visibilities and then pass the values
# and indices to a GCF_holder

def create_gcf_holder(uu, vv, vis):
    """Return GcfHolder object for a given dataset of (u,v) visibilities

    Parameters
    __________
    uu: 1D array of u coordinates from data
    vv: 1D array of v coordinates from data
    vis: ModelVisibility object
    """
    nu = uu.shape[0]
    nv = vv.shape[0]
    
    # recall that in the vis, we are padded by 2 zeros on either side in both u and v
    npix_u = vis.uu.shape[0]
    npix_v = vis.vv.shape[0]
    dense_grid_gcf = calc_dense_grid_gcf()

    # 1. Find the nearest pixel points in the fft'd image, and push to the index array
    # special case of a 1d problem means no delta u/v:
    iu0 = np.zeros(nu).astype(int)+2
    iv0 = np.zeros(nv).astype(int)+2

    if npix_u > 5:
        du = vis.uu[3] - vis.uu[2]
        iu0 = npix_u/2 + np.round((uu/du)).astype(int)

    if npix_v > 5:
        dv = vis.vv[3] - vis.vv[2]
        iv0 = npix_v/2 + np.round((vv/dv)).astype(int)

    index_arr = np.transpose(np.array((iu0, iv0)))

    # 2. Find the relative distance to this point (should be -0.5 du/v < val < 0.5 du/v) 
    u0 = uu - vis.uu[iu0]
    v0 = vv - vis.vv[iv0]

    # 3. Now we take the relative distance and find the nearest dense grid point
    iu0_grid = np.ones(nu).astype(int)*500
    iv0_grid = np.ones(nv).astype(int)*500    

    if npix_u > 5:
        iu0_grid = 500-(1001*u0/du).astype(int)

    if npix_v > 5:
        iv0_grid = 500-(1001*v0/dv).astype(int)

    # 4. Pull the gcf vals for the nearest 5 pixels around this dense grid point
    uw = dense_grid_gcf[iu0_grid,:]
    vw = dense_grid_gcf[iv0_grid,:]

    if npix_u == 5:
        uw = np.array([[0,0,1,0,0]]).astype(float)

    if npix_v == 5:
        vw = np.array([[0,0,1,0,0]]).astype(float)

    # 5. Calculate the weights (outer product) and an array of the total weights 
    all_weights = np.einsum('...a,...b->...ab', vw, uw)
    w_arr = np.einsum('...ab->...', all_weights)

    return GcfHolder(index_arr, all_weights, w_arr, uu, vv)




# The interpolation call will first calculate the gcf holder (if not provided), and then use
# those values to calculate the interpolated visibilities. 

def interpolate_uv(uu, vv, vis, gcf_holder=None):
    """Calculate interpolated visibilities

    Parameters
    __________
    uu: 1D array of u coordinates from data
    vv: 1D array of v coordinates from data
    vis: ModelVisibility object
    gcf_holder: gcf_holder object. If None, gcf_holder object will be computed. 

    Returns
    _______
    interp_vis: 2D array of interpolated visibilities with shape (N visibilities, M channels)
    gcf_holder: gcf_holder object (possibly useful for future interpolations where the gcf calculation represents a significant fraction of the computation time).
    """
    # create gcf_holder if one isn't provided
    if not gcf_holder:
        gcf_holder = create_gcf_holder(uu, vv, vis)

    # create the new interpolated visibility holder
    nvis = vv.shape[0]
    interp_vis = np.zeros((nvis, vis.freqs.shape[0]), dtype='complex')
    
    # read out the pixel indices for windowing
    u0 = gcf_holder.index_arr[:, 0].astype(int)-2
    v0 = gcf_holder.index_arr[:, 1].astype(int)-2

    # iterate through the channels and multiply weights by the windowed pixels
    # note that we could theoretically vectorize by channel as well, but this
    # produces a memory error - maybe the array is too large?
    for f in range(vis.freqs.shape[0]):
        VV_chan = vis.VV[:,:,f]
        VV = as_strided(VV_chan, shape=(VV_chan.shape[0]-4, VV_chan.shape[1]-4, 5, 5), strides=VV_chan.strides * 2)
        interp_vis[:,f] = np.einsum("ijk,ijk->i", VV[v0,u0], gcf_holder.gcf_arr)/gcf_holder.w_arr

    return interp_vis, gcf_holder
