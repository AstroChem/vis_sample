# Here we define all of the necessary functions to FFT the model image and phase shift visibilities
import numpy as np
from constants import *
from numpy.fft import fftfreq, rfftfreq, fftshift, fft2
from gridding import *
from classes import *


def transform(img):
    """Transform SkyImage object using FFT and return ModelVisibility object"""

     # we pad to prevent edge effects during interpolation
    full_fft = np.zeros(((img.data.shape[0]+4),(img.data.shape[1]+4), img.data.shape[2]), dtype='complex')

    # convert ra and dec in [arcsec] to radians, and then take the sin to
    # convert to ll, mm

    ll = np.sin(img.ra*arcsec)
    mm = np.sin(img.dec*arcsec)

    # number of elements in each array
    nl = ll.shape[0]
    nm = mm.shape[0]

    # find the spacing between the elements, as long as we're not 1D
    if nl > 1:
        dl = abs(ll[1] - ll[0]) # [radians]
    else:
        dl = 1.0

    if nm > 1:
        dm = abs(mm[1] - mm[0]) # [radians]
    else:
        dm = 1.0

    # determine uv plane coordinates in lambda
    uu = fftshift(fftfreq(nl, dl))
    vv = fftshift(fftfreq(nm, dm))

    for chan in range(img.freqs.shape[0]):
        # properly pack the data for input using fftshift to move the component at
        # RA=0,DEC=0 to the first array element: data[0,0]
        chan_fft = fftshift(fft2(fftshift(img.data[:, :, chan])))
        full_fft[2:-2,2:-2,chan] = chan_fft

        # we pad to prevent edge effects during interpolation
        uu_pad = np.lib.pad(uu, (2,2), 'constant')
        vv_pad = np.lib.pad(vv, (2,2), 'constant')

    return ModelVisibility(full_fft, uu_pad, vv_pad, img.freqs)



def phase_shift(vis, uu, vv, mu_RA, mu_DEC):
    """Shift visibilities given an offset from the phase center

    Parameters
    __________
    vis: 2D array of visibilities with shape (n visibilities, m channels)
    uu: 1D array of u coordinates from data
    vv: 1D array of v coordinates from data
    mu_RA: Offset of right ascension from phase center (arcsec)
    mu_DEC: Offset of right ascension from phase center (arcsec)

    Returns
    _______ 
    vis: 2D array of phase-shifted visibilities with shape (n visibilities, m channels)
    """
    # calculate the phase shift for each visibility
    shifts = np.exp(-2*pi*1.0j * (uu*mu_RA*arcsec + vv*mu_DEC*arcsec))

    # Go through each visibility and apply the phase shift
    vis = vis * shifts[:,np.newaxis]

    return vis


