import numpy as np
import sys
from gridding import *
from classes import *
from transforms import *
from interpolation import interpolate_uv
from file_handling import *
import time

def vis_sample(imagefile=None, uvfile=None, uu=None, vv=None, mu_RA=0, mu_DEC=0, src_distance=None, gcf_holder=None, corr_cache=None, mode="interpolate", outfile=None, verbose=False, return_gcf=False, return_corr_cache=False):
    """Sample visibilities from a sky-brightness image

    vis_sample allows you to sample visibilities from a user-supplied sky-brightness image. 

    (u,v) grid points can either be supplied by the user, or can be retrieved from a template uvfits file / measurement set.

    The results can be output either to a uvfits file or returned back to the user (for scripting)


    Parameters
    __________
    imagefile : the input sky brightness image, needs to be in a valid FITS format with units of DEG for the RA and DEC, a RADMC3D image.out file (ascii format), or a SkyImage object (use with caution)

    for uv points use:
        uvfile - uvfits file or measurement set with visibilities that the sky brightness will be interpolated to
      OR        
        uu, vv - numpy arrays - they need to be in units of lambda (i.e. number of wavelengths)

    mu_RA - (optional, default = 0) right ascension offset from phase center in arcseconds (i.e. visibilities are sampled as if the image is centered at (mu_RA, mu_DEC)

    mu_DEC - (optional, default = 0) declination offset from phase center in arcseconds (i.e. visibilities are sampled as if the image is centered at (mu_RA, mu_DEC)

    src_distance - distance to source in parsecs - only required for RADMC3D input images
 
    gcf_holder - (optional) gcf_holder object returned by previous call to vis_sample (see below return_gcf). 
        If you use this option DO NOT feed in a uvfile or uu, vv arrays, this option is intended for batch processing

    corr_cache - (optional) 2D corr_cache array output by previous call to vis_sample (see below return_corr_cache).

    mode - (optional, default = interpolate) switch between returning the interpolated output ("interpolate") or the residuals when subtracted from the data ("diff")
        NOTE - all other options are retained, so if using outfile and diff, a residual measurement set will be output
        diff mode requires an input dataset, so it will NOT work with only a gcf_holder input, and the number of channels must match
        diff mode is NOT recommended for chi sq calculation in an MCMC code, as the data will be read in for each likelihood function call (slow)
        A better solution is to use the import functions in file_handling.py to read in the data once, then use the ouput of interpolate mode to calculate chi sq

    outfile - (optional) name of output file, needs to have either a .uvfits or .ms extension consistent with extension of uvfile

    verbose - (boolean) flag to print all progress output and timing

    return_gcf - (boolean) flag to return the gcf cache to allow faster interpolation for many models   
    
    return_corr_cache - (boolean) flag to return the correction function cache to allow faster interpolation for many models 


    Usage::

    First, import the vis_sample command 
    >> from vis_sampler import vis_sample 


    1. sample my_model.fits using data.ms (u,v) points and output to interp.ms
    >> vis_sample(imagefile="my_model.fits", uvfile="data.ms", outfile='interp.ms') 


    2. sample my_model.fits using data.ms (u,v) points, interp_vis stores visibilities
    >> interp_vis = vis_sample(imagefile="my_model.fits", uvfile="data.ms") 


    In the second usage, interp_vis is the "raw" output visibility, ie just a numpy array of size [n_visibilities, n_chans]
    We can also output the caches for faster future usage:


    3. a) sample my_model.fits using data.ms (u,v) points, also store the gcf_holder
    >> interp, gcf_holder = vis_sample(imagefile="my_model.fits", uvfile="data.ms", return_gcf = True)

    3. b) sample a second model using the same (u,v) points, this is faster now
    >> interp2 = vis_sample(imagefile="second_model.fits", gcf_holder = gcf_holder)                                 


    4. a) sample my_model.fits using data.ms (u,v) points, also store the gcf_holder and corr_cache
    >> interp, gcf_holder, corr_cache = vis_sample(imagefile="my_model.fits", uvfile="data.uvfits", return_gcf = True, return_corr_cache = True)

    4. b) sample a second model using the same (u,v) points, this is faster now
    >> interp2 = vis_sample(imagefile="second_model.fits", gcf_holder = gcf_holder, corr_cache=corr_cache)
    """

    # Error cases #
    if not imagefile:
        print "Please supply an input imagefile to FFT and sample"
        return 

    if outfile:
        if return_gcf:
            print "Can only return gcf cache when not writing out to a file"
            return 
        if return_corr_cache:
            print "Can only return corr cache when not writing out to a file"
            return 
        if not uvfile: 
            print "Can only write out when there is an input data file (to clone for header info)"
            return 

    if (gcf_holder and mode=="diff"):
        print "diff mode only valid when a uvfile is supplied (to calculate residuals)"
        return 


    ###########################
    #   uu and vv retrieval   #
    ###########################

    # if we don't have uu and vv specified, then read from the cache being fed in
    if gcf_holder:
        gcf_holder = gcf_holder

    # or read them in from the data file
    elif uvfile:
        if verbose:
            print "Reading data file to interpolate onto: "+uvfile
            t0 = time.time()

        try:
            data_vis = import_data_uvfits(uvfile)
        except IOError:
            try:
                data_vis = import_data_ms(uvfile)
            except RuntimeError:
                print "Not a valid data file for interpolation. Please check that the file is a uvfits file or measurement set"
                sys.exit(1)

        if verbose: 
            t1 = time.time()
            print "Read data file to interpolate onto: "+uvfile
            print "Data read time = " + str(t1-t0)

    # if we didn't catch anything, something went wrong
    elif uu is None or vv is None:
        print "Please supply either a uvfits file to interpolate onto (uvfile), a list of uv points (uu & vv), or a gcf_holder cache"
        return

        


    ######################
    #   Read the model   #
    ######################

    # now that we have either a data_vis or a list of uu,vv points, let's import the model file
    if verbose:
        print "Reading model file: "+imagefile
        t0 = time.time()

    if isinstance(imagefile, SkyImage):
        mod_sky_img = imagefile
    elif "image.out" in imagefile:
        if src_distance is None:
             print "A source distance in pc needs to be provided in order to process a RADMC3D image file"
             return 
        else: mod_sky_img = import_model_radmc(src_distance, imagefile)
    elif "fits" in imagefile:
        mod_sky_img = import_model_fits(imagefile)
    else:
        print "Not a valid model image option. Please provide a FITS file, a RADMC3D image file, or a SkyImage object)."
        return 

    # since we clone the data file for write-out, the number of channels need to match the model
    if uvfile and (len(mod_sky_img.freqs)!=len(data_vis.freqs)):
        if mode=="diff":
            print "Number of channels in data does not match number of channels in model image. diff mode cannot continue."
            return
        else:
            print "WARNING: Number of channels in data does not match number of channels in model image. Interpolation can be completed, but model visibilities may not be able to be written to file."




    #################################
    #   Apply correction function   #
    #################################

    # necessary to correct for the effects of the convolution kernel used for interpolation
    if verbose:
        t1 = time.time()
        print "Read model file to be interpolated: "+imagefile
        print "Model read time = " + str(t1-t0)
        print "Applying corrfun"
        t0 = time.time()

    corr_cache = apply_corrfun(mod_sky_img, corr_cache=corr_cache)

    if verbose: 
        t1 = time.time()
        print "corr_fun apply time = " + str(t1-t0)




    #####################
    #   FFT the image   #
    #####################

    if verbose:
        print "Starting model FFT"
        t0 = time.time()

    mod_fft = transform(mod_sky_img)

    if verbose: 
        t1 = time.time()
        print "Model FFT complete"
        print "Model FFT time = " + str(t1-t0)
        print "Starting interpolation"



    ###################################
    #   Do the actual interpolation   #
    ###################################

    t0 = time.time()

    # dummy exists because interpolate_uv always returns a gcf_holder, but we don't need it
    if gcf_holder:
        interp, dummy = interpolate_uv(gcf_holder.uu, gcf_holder.vv, mod_fft, gcf_holder=gcf_holder)

    elif uvfile:
        interp, gcf_holder = interpolate_uv(data_vis.uu, data_vis.vv, mod_fft)

    else:
        interp, gcf_holder = interpolate_uv(uu, vv, mod_fft)

    t1 = time.time()
    if verbose: 
        print "Interpolation complete"
        print "interpolation time = " + str(t1-t0)



    ###################
    #   Phase shift   #
    ###################

    if verbose: 
        print "Starting phase shift"
        t0 = time.time()

    # calculate the pixel conversion factors - maybe this should be stored in the image class instead?         
    delt_ra = np.abs(mod_sky_img.ra[1] - mod_sky_img.ra[0])
    delt_dec = np.abs(mod_sky_img.dec[1] - mod_sky_img.dec[0])

    # we need to correct for the input image being symmetric, as the DFT expects a half pixel center offset
    # see http://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.fftfreq.html
    # additionally shift by mu_RA and mu_DEC
    interp = phase_shift(interp, gcf_holder.uu, gcf_holder.vv, mu_RA + 0.5*delt_ra, mu_DEC + 0.5*delt_dec)
        
    if verbose: 
        t1 = time.time()
        print "Phase shift complete"
        print "Phase shift time = " + str(t1-t0)



    ###################
    #   Conjugation   #
    ###################

    # conjugate the interpolated visibility to make it compatible with ALMA/SMA data
    # this has to do with an (A, B) vs (B, A) antenna/baseline convention
    # TODO - add in option to allow for other conventions?
    interp = np.conj(interp)



    #################
    #   Residuals   #
    #################

    # if in diff mode, calculate the residuals (still called interp so that the output code is the same either way)
    if mode=="diff":
        interp = data_vis.VV - interp



    ########################
    #   Now return stuff   #
    ########################

    # simplest case is just writing to a file:
    if outfile:
        # check if our ms file matches the length of the model
        if len(mod_sky_img.freqs)==len(data_vis.freqs):
            if verbose:
                print "Writing out to file: "+outfile
            interp_vis = Visibility(interp, data_vis.uu, data_vis.vv, np.ones(interp.shape), data_vis.freqs)

            # check to see what type of file we're cloning and exporting
            if "fits" in outfile:
                export_uvfits_from_clone(interp_vis, outfile, uvfile)
            elif "ms" in outfile:
                export_ms_from_clone(interp_vis, outfile, uvfile)

            # and we're done!
            return 

        # if not, then we'll center and zero-pad the model with a warning
        elif len(mod_sky_img.freqs) < len(data_vis.freqs):
            print "WARNING: number of model channels did not match number of channels in ms file. Still writing out to file, centering and zero-padding model: "+outfile
            
            # pad interp array
            npad = len(data_vis.freqs) - len(mod_sky_img.freqs)
            interp_pad = np.pad(interp, ((0,0), (int(np.floor(npad/2.)), int(np.ceil(npad/2.)))), "constant", constant_values = 0.+0.j)

            interp_vis = Visibility(interp_pad, data_vis.uu, data_vis.vv, np.ones(interp_pad.shape), data_vis.freqs)

            # check to see what type of file we're cloning and exporting
            if "fits" in outfile:
                export_uvfits_from_clone(interp_vis, outfile, uvfile)
            elif "ms" in outfile:
                export_ms_from_clone(interp_vis, outfile, uvfile)

            # and we're done!
            return
    
        else:
            print "WARNING: number of model channels greater than number of channels in ms file. Not able to write to file, returning the interpolated data."
            

    # otherwise we're going to return the raw output of the interpolation, possibly with the caches
    if return_gcf:
        if return_corr_cache:
            return interp, gcf_holder, corr_cache
        else:
            return interp, gcf_holder
    elif return_corr_cache:
        return interp, corr_cache
    else:
        return interp

