import numpy as np
import astropy.io.fits as pyfits
from constants import *
from transforms import *

# imports data from uvfits file and exports a data visibility
def import_data_uvfits(filename):
    dat = pyfits.open(filename)
    data = dat[0].data
    dhd = dat[0].header


    freq_start = dhd['CRVAL4']
    mid_chan_freq = dhd['CRPIX4']
    delt_freq = dhd['CDELT4']

    restfreq = freq_start + (mid_chan_freq-1)*delt_freq

    data_VV_raw = np.squeeze(data['data'])
    data_uu = np.squeeze(data['UU'])*restfreq
    data_vv = np.squeeze(data['VV'])*restfreq

    data_real = (data_VV_raw[:,:,0,0] + data_VV_raw[:,:,1,0])/2.
    data_imag = (data_VV_raw[:,:,0,1] + data_VV_raw[:,:,1,1])/2.
    data_wgts = (data_VV_raw[:,:,0,2] + data_VV_raw[:,:,1,2])/2.

    data_VV = data_real+data_imag*1.0j

    return Visibility(data_VV, data_uu, data_vv, data_wgts, np.arange(data_VV.shape[1])*dat[1].data['ch width'][0]/1e6), dhd


def import_model(filename):
    mod = pyfits.open(filename)
    mod_data = np.rollaxis(mod[0].data, 0, 3)
    mhd = mod[0].header

    npix_ra = mhd['NAXIS1']
    mid_pix_ra = mhd['CRPIX1']
    delt_ra = mhd['CDELT1']

    npix_dec = mhd['NAXIS2']
    mid_pix_dec = mhd['CRPIX2']
    delt_dec = mhd['CDELT2']

    nchan_vel = mhd['NAXIS3']
    mid_chan_vel = mhd['CRPIX3']
    delt_vel = mhd['CDELT3']

    # the assumption is that the RA and DEC are given in degrees, convert to arcsec 
    mod_ra = (np.arange(npix_ra)-(mid_pix_ra-1))*delt_ra*3600
    mod_dec = (np.arange(npix_dec)-(mid_pix_dec-1))*delt_dec*3600
    mod_vels = (np.arange(nchan_vel)-(mid_chan_vel-1))*delt_vel

    return SkyImage(mod_data, mod_ra, mod_dec, mod_vels)

# ONLY CAN CLONE UVFITS
# TODO - FIGURE OUT HOW TO WRITE FROM SCRATCH
def export_uvfits_from_clone(vis, outfile, uvfits_clone):
    clone = pyfits.open(uvfits_clone)
    clone_data = clone[0].data

    data_array = np.zeros([vis.VV.shape[0], vis.VV.shape[1], 2, 3])
    data_array[:,:,0,0] = np.real(vis.VV)
    data_array[:,:,1,0] = np.real(vis.VV)
    data_array[:,:,0,1] = np.imag(vis.VV)
    data_array[:,:,1,1] = np.imag(vis.VV)
    data_array[:,:,0,2] = vis.wgts
    data_array[:,:,1,2] = vis.wgts

    clone_data['data'] = np.expand_dims(np.expand_dims(np.expand_dims(data_array, 1),1),1)

    clone.writeto(outfile)
    print "Wrote " + outfile
