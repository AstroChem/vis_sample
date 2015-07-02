import numpy as np
from numpy.fft import fftfreq, rfftfreq, fftshift, fft2
from gridding import *
from classes import *
from transforms import *
import time
from scipy import weave
from scipy.weave import converters


def calc_dense_grid_gcf():
    dense_grid = np.linspace(-0.5, 0.5, 1000)
    dense_grid_gcf = np.zeros((1000,6))
    for i, grid_pnt in enumerate(dense_grid):
        eta = (np.arange(-3,3) + grid_pnt)/3.0
        dense_grid_gcf[i] = gcffun(eta, 1.0)
    return dense_grid_gcf

def interpolate_uv_single_point_dense_grid(u, v, l, vis, dense_grid_gcf):
    # Note that vis.uu goes from positive to negative (East-West)
    # and vis.vv goes from negative to positive (North-South)

    #t00 = time.time()
    #t0 = time.time()
    # 1. Find the nearest gridpoint in the FFT'd image.
    uvcnt = vis.uu.shape[0]
    # determine the uu and vv distance for grid points
    du = vis.uu[1] - vis.uu[0]
    dv = vis.vv[1] - vis.vv[0]

    iu0 = uvcnt/2+int(u/du+0.5)
    iv0 = uvcnt/2+int(v/dv+0.5)

    # now find the relative distance to this nearest grid point (not absolute)
    u0 = u - vis.uu[iu0]
    v0 = v - vis.vv[iv0]

    # 2. Calculate the appropriate u and v grid points
    # find the nearest dense grid point
    iu0_grid = 500+int(1000*u0/du)
    iv0_grid = 500+int(1000*v0/dv)

    VV = vis.VV[(iv0-3):(iv0+3), (iu0-3):(iu0+3), l] # Array is packed like the image

    #t1 = time.time()
    #print "total search time = " + str(t1-t0)

    # 3. Calculate the weights corresponding to these 6 nearest pixels (gcffun)
    #t0 = time.time()
    uw = dense_grid_gcf[iu0_grid,:]
    vw = dense_grid_gcf[iv0_grid,:]
    #t1 = time.time()
    #print "gcf calc time = " + str(t1-t0)

    all_weights = np.outer(uw, vw)

    # 4. Normalization such that it has an area of 1. Divide by w later.
    w = (uw[0]+uw[1]+uw[2]+uw[3]+uw[4]+uw[5])*(vw[0]+vw[1]+vw[2]+vw[3]+vw[4]+vw[5])

    # 5. Loop over all 36 grid indices and sum to find the interpolation.
    cum = np.sum(all_weights*VV)/w

    #t11 = time.time()
    #print "total single point time = " + str(t11-t00)

    return cum, np.array([(iu0-3), (iu0+3), (iv0-3), (iv0+3)]), all_weights, w



# called ModGrid in gridding.c (KR code) and in Model.for (MIRIAD)
# Uses spheroidal wave functions to interpolate a model to a (u,v) coordinate.
# u,v are in [klam]
def interpolate_uv_single_point(u, v, l, vis):
    # Note that vis.uu goes from positive to negative (East-West)
    # and vis.vv goes from negative to positive (North-South)

    #t00 = time.time()
    #t0 = time.time()
    # 1. Find the nearest gridpoint in the FFT'd image.
    iu0 = np.argmin(abs(u - vis.uu))
    iv0 = np.argmin(abs(v - vis.vv))
    #t1 = time.time()
    #print "nearest neighbor search time = " + str(t1-t0)

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
    #t0 = time.time()
    uw = gcffun(etau, 1.0)
    vw = gcffun(etav, 1.0)
    all_weights = np.outer(vw, uw)
    #t1 = time.time()
    #print "gcf calc time = " + str(t1-t0)    

    # 4. Normalization such that it has an area of 1. Divide by w later.
    w = np.sum(all_weights)

    # 5. Loop over all 36 grid indices and sum to find the interpolation.
    cum = 0.+0.j
    for i in range(6):
        for j in range(6):
            cum += uw[i] * vw[j] * VV[j,i,l] # Array is packed like the image

    cum = cum/w
    
    #time11 = time.time()
    #print "total single point time = " + str(t11-t00)

    return cum, np.array([uind[0], (uind[5]+1), vind[0], (vind[5]+1)]), all_weights, w


# interpolate a model for many u,v points and store the weights for later
def interpolate_uv_initial(uu, vv, vis):
    uvcnt = vv.shape[0]
    interp_vis = np.zeros((uvcnt), dtype='complex')
    index_arr = np.zeros((uvcnt, 4))
    gcf_arr = np.zeros((uvcnt, 6, 6))
    w_arr = np.zeros((uvcnt))

    dense_grid_gcf = calc_dense_grid_gcf()

    print "interpolating " +str(uvcnt) + " visibilities"
    for i in range(uvcnt):
#        interp_vis[i], index_arr[i], gcf_arr[i], w_arr[i] = interpolate_uv_single_point(uu[i], vv[i], 0, vis)
         interp_vis[i], index_arr[i], gcf_arr[i], w_arr[i] = interpolate_uv_single_point_dense_grid(uu[i], vv[i], 0, vis, dense_grid_gcf)

    return interp_vis, GCF_holder(index_arr, gcf_arr, w_arr, uu, vv)
    
# interpolate a model for many u,v points
def interpolate_uv(uu, vv, vis, gcf_holder='', return_gcf=False):
    t0 = time.time()
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
                interp_vis[i,l] = np.sum(gcf_holder.gcf_arr[i]*VV)/gcf_holder.w_arr[i]
        if (return_gcf==True):
            t1 = time.time()
            print "total interpolation time = " + str(t1-t0)
            return interp_vis, gcf_holder
        else:
            t1 = time.time()
            print "total interpolation time = " + str(t1-t0)
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
                interp_vis[i,l] = np.sum(gcf_holder.gcf_arr[i]*VV)/gcf_holder.w_arr[i]
        t1 = time.time()
        print "total interpolation time = " + str(t1-t0)
        return interp_vis
