import numpy as np
from constants import *
from numpy.fft import fftfreq, rfftfreq, fftshift, fft2
from gridding import *
from classes import *

# Given a new model centroid in the image plane (in arcseconds), shift the
# visibilities by corresponding amount
def phase_shift(vis, mu_RA, mu_DEC):
    mu = np.array([mu_RA, mu_DEC]) * arcsec # [radians]
    nu = vis.uu.shape[0]
    nv = vis.vv.shape[0]
    
    # Go through each visibility and apply the phase shift
    for l in range(vis.lams.shape[0]):
        for i in range(nu):
            for j in range(nv):
                # Convert from [klam] to [lam]
                R = np.array([vis.uu[i], vis.vv[j]]) * 1e3 #[lam]
                # Not actually in polar phase form
                shift = np.exp(-2*pi*1.0j * np.sum(R*mu))
                vis.VV[j,i,l] = vis.VV[j,i,l] * shift



# Transform the SkyImage using FFT
def transform(img):
    full_out = np.zeros(img.data.shape, dtype='complex')

    # convert ra and dec in [arcsec] to radians, and then take the sin to
    # convert to ll, mm
    ll = np.sin(img.ra * arcsec)
    mm = np.sin(img.dec * arcsec)

    # number of elements in each array
    nl = ll.shape[0]
    nm = mm.shape[0]

    # find the spacing between the elements
    dl = abs(ll[1] - ll[0]) # [radians]
    dm = abs(mm[1] - mm[0]) # [radians]

    # determine uv plane coordinates in klam
    uu = fftshift(fftfreq(nl, dl)) * 1e-3 # [klam]
    vv = fftshift(fftfreq(nm, dm)) * 1e-3 # [klam]

    for index in range(img.lams.shape[0]):
        data = img.data[::-1, :, index]

        # properly pack the data for input using fftshift to move the component at
        # RA=0,DEC=0 to the first array element: data[0,0]

        chan_out = fftshift(fft2(fftshift(data)))
        full_out[:,:,index] = chan_out

    return Model_Visibility(full_out, uu, vv, img.lams)
