# Here we calculate the gridding convolution function (gcf) and the correct function (corrfun)
# for a given distance in grid space (eta). We interpolate over the nearest 7 pixels.
# The gcf is used for the interpolation later, and the corrfun is applied to the image prior
# to the FFT to correct for the future convolution of the gcf (ie corrfun is the FT'd inverse)

# TODO: correct etalim for using 7 pixels rather than 6

import numpy as np
from scipy import weave
from scipy.weave import converters
import time



# First we definte the prolate spheroidal functions used for calculating the gcf and corrfun.
# This funtion definition comes from Schwabb's derivations

def spheroid_weave(eta=0):
    etalim = 0.75
    eta = abs(eta)

    if eta <= etalim:
        nn = eta**2 - etalim**2
        return spheroid_weave_1(nn)
    elif eta <= 1.00000000001:
        nn = eta**2 - 1.0
        return spheroid_weave_2(nn)

    return 1e30


def spheroid_weave_1(nn):
    support = """
        #include <math.h>
        #include <stdlib.h>
        """
    code = """
        double n = (double) nn;
        return_val = (pow(0.2312756*n,4) + pow(-0.5335581*n,3) + pow(0.627866*n,2) + -0.3644705*n + 0.08203343)/(pow(0.2078043*n,2) + 0.8212018*n + 1.0);
    """
    return weave.inline(code, arg_names=['nn'], support_code = support, libraries = ['m'], type_converters=converters.blitz)


def spheroid_weave_2(nn):
    support = """
        #include <math.h>
        #include <stdlib.h>
        """
    code = """
        double n = (double) nn;
        return_val = (pow(0.06412774*n,4) + pow(-0.1201436*n,3) + pow(0.1021332*n,2) + -0.03697768*n + 0.004028559)/(pow(0.2918724*n,2) + 0.9599102*n + 1.0);
    """
    return weave.inline(code, arg_names=['nn'], support_code = support, libraries = ['m'], type_converters=converters.blitz)




# now we define the functions that will calculate the gcf
def gcf_single(eta):
    return ((abs(1 - eta**2)))*spheroid_weave(eta)

def gcffun(etas):
    return [gcf_single(eta) for eta in etas]





# Finally we need a function to apply the corrfun to a sky brightness image
#
# Caching the corrfun values results in a ~40% speed increase for future model corrections
# corr_cache can be returned through the return_cache flag, then fed back in later

def apply_corrfun(img, mu_ra, mu_dec, corr_cache="", return_cache=False):
    ndec, nra, nvel = img.data.shape
    if (corr_cache==""):
        del_ra = abs(img.ra[1] - img.ra[0])
        del_dec = abs(img.dec[1] - img.dec[0])
        maxra = del_ra * nra/2
        maxdec = del_dec * ndec/2

        # If the image will be later offset via a phase shift, then this means that
        # the corrfun will need to be applied *as if the image were already offset*
        corr_cache = np.zeros([ndec, nra])
        eta_x = (img.ra + mu_ra)/maxra
        # not sure why -del_dec was applied here before - has to do with pixel shift I think
        eta_y = (img.dec -del_dec + mu_dec)/maxdec

        spheroid_vectorized = np.vectorize(spheroid_weave)
        corr_x = 1.0/spheroid_vectorized(eta_x)
        corr_y = 1.0/spheroid_vectorized(eta_y)

        corr_cache = np.outer(corr_x, corr_y)


    for k in range(nvel):
        img.data[:,:,k] = img.data[:,:,k]*corr_cache

    if (return_cache==True):
        return corr_cache




