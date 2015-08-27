# Here we calculate the gridding convolution function (gcf) and the correct function (corrfun)
# for a given distance in grid space (eta). We interpolate over the nearest 5 pixels.
# The gcf is used for the interpolation later, and the corrfun is applied to the image prior
# to the FFT to correct for the future convolution of the gcf (ie corrfun is the FT'd inverse)

import numpy as np
from scipy import weave
from scipy.weave import converters
import time



# First we definte the prolate spheroidal functions used for calculating the gcf and corrfun.
# This funtion definition comes from Schwab's derivations

def spheroid_weave(eta=0):
    """
    Calculate value of spheroidal function used for gridding convolution function at a given value eta

    The rational approximation for this function is presented in F.R. Schwab's "Optimal Gridding of Visibility Data in Radio Interferometry" in Indirect Imaging: Measurement and Processing for Indirect Imaging, 1984
    m=5, alpha = 1

    Parameter
    _________
    eta: float between -1 and 1 

    Returns
    _______
    out: float. Value of spheroidal function at eta
    
    """
    nn = eta**2 - 1**2

    if np.abs(eta) < 1.0000000000001:
        support = """
            #include <math.h>
            #include <stdlib.h>
            """
        code = """
            double n = (double) nn;
            return_val = (0.01624782*pow(n,6) + -0.05350728*pow(n,5) + 0.1464354*pow(n,4) + -0.2347118*pow(n,3) + 0.2180684*pow(n,2) + -0.09858686*n + 0.01466325)/(0.2177793*n + 1);
        """
        return weave.inline(code, arg_names=['nn'], support_code = support, libraries = ['m'], type_converters=converters.blitz)

    else:
        return 1e30


# now we define the functions that will calculate the gcf
def gcf_single(eta):
    return (abs(1 - eta**2))*spheroid_weave(eta)

def gcffun(etas):
    return [gcf_single(eta) for eta in etas]





# Finally we need a function to apply the corrfun to a sky brightness image
#
# Caching the corrfun values results in a ~40% speed increase for future model corrections
# corr_cache can be returned through the return_cache flag, then fed back in later

def apply_corrfun(img, corr_cache=None):
    """Applies correction function to sky brightness image

    Parameters
    __________
    img: SkyImage object (required)
    corr_cache: ndarray (optional) 2D array that is the same shape as img.data (i.e., the number of elements in each dimension is equal to the number of pixels along that dimension in the original image file). Feeding this array in as a parameter avoids having to recalculate the correction function in vis_sample

    Returns
    _______
    corr_cache: ndarray. 2D array that is the same shape as img.data (i.e., the number of elements in each dimension is equal to the number of pixels along that dimension in the original image file). Feeding this array in as a parameter avoids having to recalculate the correction function in vis_sample
    
    """

    ndec, nra, nvel = img.data.shape

    if (corr_cache is None):
        corr_cache = np.zeros([ndec, nra])

        eta_x = np.array([0])
        eta_y = np.array([0])

        if nra > 1:
            del_ra = img.ra[1] - img.ra[0]
            maxra = abs(del_ra) * nra/2
            eta_x = (img.ra)/maxra

        if ndec > 1:
            del_dec = img.dec[1] - img.dec[0]
            maxdec = abs(del_dec) * ndec/2
            eta_y = (img.dec)/maxdec

        spheroid_vectorized = np.vectorize(spheroid_weave)
        corr_x = 1.0/spheroid_vectorized(eta_x)
        corr_y = 1.0/spheroid_vectorized(eta_y)

        corr_cache = np.outer(corr_y, corr_x)

    for k in range(nvel):
        img.data[:,:,k] = img.data[:,:,k]*corr_cache

    return corr_cache




