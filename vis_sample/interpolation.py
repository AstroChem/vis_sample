import numpy as np
from numpy.fft import fftfreq, rfftfreq, fftshift, fft2
from gridding import *
from classes import *
from transforms import *

# called ModGrid in gridding.c (KR code) and in Model.for (MIRIAD)
# Uses spheroidal wave functions to interpolate a model to a (u,v) coordinate.
# u,v are in [klam]
def interpolate_uv_single_point(u, v, l, vis):
    # Note that vis.uu goes from positive to negative (East-West)
    # and vis.vv goes from negative to positive (North-South)

    # 1. Find the nearest gridpoint in the FFT'd image.
    iu0 = np.argmin(abs(u - vis.uu))
    iv0 = np.argmin(abs(v - vis.vv))

    # now find the relative distance to this nearest grid point (not absolute)
    u0 = u - vis.uu[iu0]
    v0 = v - vis.vv[iv0]

    # determine the uu and vv distance for 3 grid points (could be later taken out)
    du = abs(vis.uu[3] - vis.uu[0])
    dv = abs(vis.vv[3] - vis.vv[0])

    # 2. Calculate the appropriate u and v indexes for the 6 nearest pixels
    # (3 on either side)
    # Are u0 and v0 to the left or the right of the index?
    # we want to index three to the left, three to the right

    if iu0 < 2 or iv0 < 2 or iu0 > (vis.uu.shape[0]-2) or iv0 > (vis.vv.shape[0]-2):
        return 0, np.array([0,6,0,6]), np.zeros(36), 1
    if u0 >= 0.0:
        # To the right of the index
        uind = np.arange(iu0-2, iu0+4)
    else:
        # To the left of the index
        uind = np.arange(iu0-3, iu0+3)

    if v0 >= 0.0:
        # To the right of the index
        vind = np.arange(iv0-2, iv0+4)
    else:
        # To the left of the index
        vind = np.arange(iv0-3, iv0+3)

    etau = (vis.uu[uind] - u)/(du)
    etav = (vis.vv[vind] - v)/(dv)

    VV = vis.VV[vind[0]:(vind[5]+1), uind[0]:(uind[5]+1)] # Array is packed like the image

    # 3. Calculate the weights corresponding to these 6 nearest pixels (gcffun)
    uw = gcffun(etau, 1.0)
    vw = gcffun(etav, 1.0)
    all_weights = np.outer(vw, uw)

    # 4. Normalization such that it has an area of 1. Divide by w later.
    w = np.sum(all_weights)

    # 5. Loop over all 36 grid indices and sum to find the interpolation.
    cum = 0.+0.j
    for i in range(6):
        for j in range(6):
            cum += uw[i] * vw[j] * VV[j,i,l] # Array is packed like the image

    cum = cum/w

    return cum, np.array([uind[0], (uind[5]+1), vind[0], (vind[5]+1)]), np.reshape(all_weights, 36), w


# interpolate a model for many u,v points and store the weights for later
def interpolate_uv_initial(uu, vv, vis):
    uvcnt = vv.shape[0]
    interp_vis = np.zeros((uvcnt), dtype='complex')
    index_arr = np.zeros((uvcnt, 4))
    gcf_arr = np.zeros((uvcnt, 36))
    w_arr = np.zeros((uvcnt))

    for i in range(uvcnt):
        interp_vis[i], index_arr[i], gcf_arr[i], w_arr[i] = interpolate_uv_single_point(uu[i], vv[i], 0, vis)
            
    return interp_vis, GCF_holder(index_arr, gcf_arr, w_arr, uu, vv)
    
# interpolate a model for many u,v points
def interpolate_uv(uu, vv, vis, gcf_holder='', return_gcf=False):
    if (gcf_holder==""):
        uvcnt = vv.shape[0]
        interp_vis = np.zeros((uvcnt, vis.lams.shape[0]), dtype='complex')
        interp_vis[:,0], gcf_holder = interpolate_uv_initial(uu, vv, vis)
        for l in range(1, vis.lams.shape[0]):
            for i in range(uvcnt):
                uind0 = gcf_holder.index_arr[i, 0]
                uind5 = gcf_holder.index_arr[i, 1]
                vind0 = gcf_holder.index_arr[i, 2]
                vind5 = gcf_holder.index_arr[i, 3]
                VV = vis.VV[vind0:vind5, uind0:uind5, l]
                interp_vis[i,l] = np.sum(gcf_holder.gcf_arr[i]*np.reshape(VV,36))/gcf_holder.w_arr[i]
        if (return_gcf==True):
            return interp_vis, gcf_holder
        else:
            return interp_vis
        
    else:
        uvcnt = vv.shape[0]
        interp_vis = np.zeros((uvcnt, vis.lams.shape[0]), dtype='complex')
        for l in range(0, vis.lams.shape[0]):
            for i in range(uvcnt):
                uind0 = gcf_holder.index_arr[i, 0]
                uind5 = gcf_holder.index_arr[i, 1]
                vind0 = gcf_holder.index_arr[i, 2]
                vind5 = gcf_holder.index_arr[i, 3]
                VV = vis.VV[vind0:vind5, uind0:uind5, l]
                interp_vis[i,l] = np.sum(gcf_holder.gcf_arr[i]*np.reshape(VV,36))/gcf_holder.w_arr[i]
        return interp_vis
