import numpy as np
from gridding import *
from classes import *
from transforms import *
from interpolation import interpolate_uv
from file_handling import *
import time

def vis_sample(imagefile=None, uvfile=None, uu=None, vv=None, mu_RA=0, mu_DEC=0, src_distance = None, gcf_holder=None, corr_cache=None, outfile=None, verbose=False, return_gcf=False, return_corr_cache=False):
    """Sample visibilities from a sky-brightness image

    vis_sample allows you to sample visibilities from a user-supplied sky-brightness image. 

    (u,v) grid points can either be supplied by the user, or can be retrieved from a template uvfits file / measurement set.

    The results can be output either to a uvfits file or returned back to the user (for scripting)

    Parameters
    __________
    imagefile : the input sky brightness image, needs to be in a valid FITS format with units of DEG for the RA and DEC or a RADMC3D image.out file (ascii format)

    for uv points use:
        uvfile - uvfits file or measurement set with visibilities that the sky brightness will be interpolated to
      OR        
        uu, vv - numpy arrays - they need to be in units of lambda (i.e. number of wavelengths)

    mu_RA - right ascension offset from phase center in arcseconds 

    mu_DEC - declination offset from phase center in arcseconds

    src_distance - distance to source in parsecs - only required for RADMC3D input images
 
    gcf_holder - (optional) gcf_holder object returned by previous call to vis_sample (see below return_gcf). 
                If you use this option DO NOT feed in a uvfile or uu, vv arrays. They will be used by default and you'll see no speed increase

    corr_cache - (optional parameter) 2D corr_cache array output by previous call to vis_sample (see below return_corr_cache). 

    outfile - name of optional output file, needs to have either a .uvfits or .ms extension consistent with extension of uvfile

    verbose - prints all progress output and timing

    return_gcf - (boolean) flag to return the gcf cache to allow faster interpolation for many models   
    
    return_corr_cache - (boolean) flag to return the correction function cache to allow faster interpolation for many models 


    Usage::

    >> from vis_sample import vis_sample                                                                            # import the vis_sample command  

    >> vis_sample(imagefile="my_model.fits", uvfile="data.uvfits", outfile='interp.uvfits')         # sample my_model using data (u,v) points and output to interp.uvfits

    >> interp_vis = vis_sample(imagefile="my_model.fits", uvfile="data.uvfits")                                     # sample my_model using data (u,v) points, interp_vis stores visibilities

    In the second usage, interp_vis is the "raw" output visibility, ie just a numpy array of size [n_visibilities, n_chans]

    We can also output the caches for faster future usage:

    >> interp, gcf_holder = vis_sample(imagefile="my_model.fits", uvfile="data.uvfits", return_gcf = True)          # sample my_model using data (u,v) points, also store the gcf_holder

    >> interp2 = vis_sample(imagefile="second_model.fits", gcf_holder = gcf_holder)                                 # sample a second model using the same (u,v) points, this is faster now

    >> interp, gcf_holder, corr_cache = vis_sample(imagefile="my_model.fits", uvfile="data.uvfits", return_gcf = True, return_corr_cache = True)          # sample my_model using data (u,v) points, also store the gcf_holder and corr_cache

    >> interp2 = vis_sample(imagefile="second_model.fits", gcf_holder = gcf_holder, corr_cache=corr_cache)                                 # sample a second model using the same (u,v) points, this is faster now
    """

    # Error cases #
    if not imagefile:
        print "Please supply an input imagefile to FFT and sample"
        return 

    if outfile:
        if return_gcf:
            print "Can only return gcf cache when not writing out to a file"
            return 
        if not uvfile: 
            print "Can only write out when there is an input data file (for header info)"
            return 

    if not (uu and vv or gcf_holder or uvfile):
         print "Please supply either a uvfits file to interpolate onto (uvfile), a list of uv points (uu & vv), or a gcf_holder cache"


    ###########################
    #   uu and vv retrieval   #
    ###########################

    # if we don't have uu and vv specified, then read them in from the data file
    if uvfile:
        if "fits" in uvfile:
            data_vis = import_data_uvfits(uvfile)
        elif "ms" in uvfile:
            data_vis = import_data_ms(uvfile)
        else:
            print "not a valid data file to interpolate onto"
        if verbose: print "Read data file to interpolate onto: "+uvfile

    # or read from the cache being fed in
    elif gcf_holder:
        gcf_holder = gcf_holder
        


    ####################################################
    #   Read the model and apply correction function   #
    ####################################################

    # now that we have either a data_vis or a list of uu,vv points, let's import the model file
    if "fits" in imagefile:
        mod_sky_img = import_model_fits(imagefile)
    elif "image.out" in imagefile:
        if src_distance is None:
             print "A source distance in pc needs to be provided in order to process a RADMC3D image file"
             return 
        else: mod_sky_img = import_model_radmc(src_distance, imagefile)

    if uvfile and (len(mod_sky_img.freqs)!=len(data_vis.freqs)):
        print "WARNING: Number of channels in data does not match number of channels in model image. Interpolation can be completed, but model visibilities cannot be written to file."

    # now apply the correction function
    if verbose:
        print "Read model file to be interpolated: "+imagefile
        t0 = time.time()
        print "Applying corrfun"

    corr_cache = apply_corrfun(mod_sky_img, corr_cache=corr_cache)

    if verbose: 
        t1 = time.time()
        print "corr_fun apply time = " + str(t1-t0)




    #####################
    #   FFT the image   #
    #####################

    if verbose:
        print "Starting FFT"
        t0 = time.time()

    mod_fft = transform(mod_sky_img)

    if verbose: 
        t1 = time.time()
        print "fft time = " + str(t1-t0)
        print "FFT complete"
        print "Starting interpolation"


    ###################################
    #   Do the actual interpolation   #
    ###################################

    t0 = time.time()

    if gcf_holder:
        interp, dummy = interpolate_uv(gcf_holder.uu, gcf_holder.vv, mod_fft, gcf_holder=gcf_holder)

    elif uvfile:
        interp, gcf_holder = interpolate_uv(data_vis.uu, data_vis.vv, mod_fft)

    else:
        interp, gcf_holder = interpolate_uv(uu, vv, mod_fft)

    t1 = time.time()
    if verbose: 
        print "interpolation time = " + str(t1-t0)
        print "Interpolation complete"


    #############################
    #   Phase shift if needed   #
    #############################

    if not (mu_RA == 0 and mu_DEC == 0):
        if verbose: 
            print "Starting phase shift"
            t0 = time.time()
        
        interp = phase_shift(interp, gcf_holder.uu, gcf_holder.vv, mu_RA, mu_DEC)
        
        if verbose: 
            t1 = time.time()
            print "Phase shift time = " + str(t1-t0)
            print "Phase shift complete"


    ########################
    #   Now return stuff   #
    ########################

    # simplest case is just writing to a file:
    if outfile and len(mod_sky_img.freqs)==len(data_vis.freqs):
        if verbose:
            print "Writing out to file: "+outfile
        interp_vis = Visibility(interp, data_vis.uu, data_vis.vv, np.ones(interp.shape), data_vis.freqs)
        if "fits" in outfile:
            export_uvfits_from_clone(interp_vis, outfile, uvfile)
        elif "ms" in outfile:
            export_ms_from_clone(interp_vis, outfile, uvfile)

        # and we're done!
        return 

    if return_gcf:
        if return_corr_cache:
            return interp, gcf_holder, corr_cache
        else:
            return interp, gcf_holder
    elif return_corr_cache:
        return interp, corr_cache
    else:
        return interp

