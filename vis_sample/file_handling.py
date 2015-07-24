import numpy as np
import astropy.io.fits as pyfits
from constants import *
from transforms import *
import sys
import shutil

def import_data_uvfits(filename):
    """Imports data from uvfits file and returns Visibility object"""
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

    return Visibility(data_VV, data_uu, data_vv, data_wgts, np.arange(data_VV.shape[0])*dat[1].data['ch width'][0]/1e6)


# CASA interfacing code comes from Peter Williams' casa-python and casa-data package
# commands for retrieving ms data are from Sean Andrews
def import_data_ms(filename):
    """Imports data from a casa measurement set (ms) and returns Visibility object"""
    try:
        import casac
    except:
        print "casac was not able to be imported, make sure all dependent packages are installed"
        print "try: conda install -c pkgw casa-python casa-data"
        sys.exit(1)

    tb = casac.casac.table()
    ms = casac.casac.ms()
    
    # Use CASA table tools to get columns of UVW, DATA, WEIGHT, etc.
    tb.open(filename)
    data    = tb.getcol("DATA")
    uvw     = tb.getcol("UVW")
    weight  = tb.getcol("WEIGHT")
    ant1    = tb.getcol("ANTENNA1")
    ant2    = tb.getcol("ANTENNA2")
    tb.close()

    # Use CASA ms tools to get the channel/spw info
    ms.open(filename)
    spw_info = ms.getspectralwindowinfo()
    start_freq = spw_info["0"]["Chan1Freq"]
    chan_width = spw_info["0"]["ChanWidth"]
    nchan = spw_info["0"]["NumChan"]
    rfreq = start_freq + nchan/2.0*chan_width
    ms.close()

    data = np.squeeze(data)

    # break out the u, v spatial frequencies, convert from m to lambda
    uu = uvw[0,:]*rfreq/(cc/100)
    vv = uvw[1,:]*rfreq/(cc/100)

    # polarization averaging
    Re_xx = data[0,:].real
    Re_yy = data[1,:].real
    Im_xx = data[0,:].imag
    Im_yy = data[1,:].imag
    weight_xx = weight[0,:]
    weight_yy = weight[1,:]

    # - weighted averages
    Re = (Re_xx*weight_xx + Re_yy*weight_yy) / (weight_xx + weight_yy)
    Im = (Im_xx*weight_xx + Im_yy*weight_yy) / (weight_xx + weight_yy)
    wgts = (weight_xx + weight_yy)

    # toss out the autocorrelation placeholders
    xc = np.where(ant1 != ant2)[0]

    data_real = Re[:,xc]
    data_imag = Im[:,xc]
    data_wgts = wgts[xc]
    data_uu = uu[xc]
    data_vv = vv[xc]

    data_VV = data_real+data_imag*1.0j

    return Visibility(data_VV, data_uu, data_vv, data_wgts, (np.arange(data_VV.shape[0])*chan_width + start_freq)/1e6)


# imports model from a FITS file - note the assumptions on dimensions
def import_model(filename):
    """Imports model from a FITS file and returns SkyImage object

    Note the assumption that RA and DEC are given in degrees (converted to arcsec when returned in SkyImage object)

    """
    mod = pyfits.open(filename)
    mod_data = np.rollaxis(mod[0].data, 0, 3)
    mhd = mod[0].header

    npix_ra = mhd['NAXIS1']
    mid_pix_ra = mhd['CRPIX1']
    delt_ra = mhd['CDELT1']
    if delt_ra < 0:
        mod_data = np.fliplr(mod_data)

    npix_dec = mhd['NAXIS2']
    mid_pix_dec = mhd['CRPIX2']
    delt_dec = mhd['CDELT2']
    if delt_dec < 0:
        mod_data = np.flipud(mod_data)

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


# ONLY CAN CLONE MS
# TODO - FIGURE OUT HOW TO WRITE FROM SCRATCH
def export_ms_from_clone(vis, outfile, ms_clone):
    try:
        import casac
    except:
        print "casac was not able to be imported, make sure all dependent packages are installed"
        print "try: conda install -c pkgw casa-python casa-data"
        sys.exit(1)

    shutil.copytree(ms_clone, outfile)

    tb = casac.casac.table()
    
    # Use CASA table tools to fill new DATA and WEIGHT
    tb.open(outfile, nomodify=False)

    # we need to pull the antennas and find where the autocorrelation values aren and aren't
    ant1    = tb.getcol("ANTENNA1")
    ant2    = tb.getcol("ANTENNA2")
    ac = np.where(ant1 == ant2)[0]
    xc = np.where(ant1 != ant2)[0]

    data_array = np.zeros((2, vis.VV.shape[1], ant1.shape[0])).astype(complex)

    # fill the xc with our interpolation
    data_array[0, :, xc] = vis.VV
    data_array[1, :, xc] = vis.VV

    # fill the ac with 0's
    data_array[0, :, ac] = 0 + 0j
    data_array[1, :, ac] = 0 + 0j

    # now do the same with the weights
    weights = np.zeros((2, ant1.shape[0]))
    weights[0, xc] = np.mean(vis.wgts, axis=1)
    weights[1, xc] = np.mean(vis.wgts, axis=1)

    tb.putcol("DATA", data_array)
    tb.putcol("WEIGHT", weights)
    tb.flush()
    tb.close()

    print "Wrote " + outfile
