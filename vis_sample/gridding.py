import numpy as np
from scipy import weave
from scipy.weave import converters
from math import *

# assumes m = 6
def corrfun_single(eta, alpha):
    return spheroid_weave(eta, alpha)

def corrfun(etas, alpha):
    return [corrfun_single(eta, alpha) for eta in etas]

def gcffun_single(eta, alpha):
    return ((abs(1 - eta**2))**alpha)*spheroid_weave(eta, alpha)

def gcffun(etas, alpha):
    return [gcffun_single(eta, alpha) for eta in etas]


# Apply the correction function to the image.
def apply_corrfun(img, alpha, mu_RA, mu_DEC, corr="", return_corr=False):
    if (corr==""):
        ndec, nra, nvel = img.data.shape
        # The size of one half-of the image.
        # sometimes ra and dec will be symmetric about 0, othertimes they won't
        # so this is a more robust way to determine image half-size
        del_ra = abs(img.ra[1] - img.ra[0])
        del_dec = abs(img.dec[1] - img.dec[0])
        maxra = del_ra * nra/2
        maxdec = del_dec * ndec/2

        # If the image will be later offset via a phase shift, then this means that
        # the corrfunction will need to be applied *as if the image were already
        # offset.*

        corr = np.zeros([ndec, nra])
        for i in range(nra):
            for j in range(ndec):
                eta_x = (img.ra[i] + mu_RA)/maxra
                # not sure why -del_dec needs to be applied - has to do with pixel shift I think
                eta_y = (img.dec[j] - del_dec + mu_DEC)/maxdec
                if (abs(eta_x) <= 1.0 and abs(eta_y) <= 1.0):
                    corr[j, i] = 1.0/(corrfun_single(eta_x, alpha) * corrfun_single(eta_y, alpha))

        for k in range(nvel):
            img.data[:,:,k] = img.data[:,:,k]*corr

        if (return_corr==True):
            return corr
    else:
        for k in range(nvel):
            img.data[:,:,k] = img.data[:,:,k]*corr


# assumes m = 6
def spheroid_weave(eta=0, alpha=1):
    # assuming alpha=1 for now

    etalim = 0.75  # specific for m = 6

    eta = abs(eta)

    if eta <= etalim:
        nn = eta**2 - etalim**2
        return spheroid_weave_1(nn)

    elif eta <= 1.00000000001:
        nn = eta**2 - 1.0
        return spheroid_weave_2(nn)

    return 0.


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

