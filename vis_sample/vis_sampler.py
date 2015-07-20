import numpy as np
from gridding import *
from classes import *
from transforms import *
from interpolation import interpolate_uv
from file_handling import *
import time


# currently the imagefile needs to be formatted with units of DEG in RA and DEC
# units for uu and vv are LAMBDA (ie number of wavelengths)
# units for mu_RA and mu_DEC are arcsec

def vis_sample(imagefile=0, uvfile=0, uu=0, vv=0, mu_RA=0, mu_DEC=0, gcf_holder=0, corr_cache=0, writefile=False, outfile="", verbose=False, return_gcf=False, return_corr_cache=False):

    # Error cases #
    if imagefile==0:
        print "Please supply an input imagefile to FFT and sample"
        return ""

    if (writefile==True and outfile==""):
        print "Please supply an output file name (outfile) if writing to a file"
        return ""

    if (writefile==True and return_gcf==True):
        print "Can only return gcf cache when not writing out to a file"
        return ""

    if (writefile==True and uvfile==0):
        print "Can only write out when there is an input data file (for header info)"
        return ""

    if (uvfile == 0):
        if (np.size(uu) == 0) or (np.size(vv) == 0):
            if (gcf_holder == 0):
                print "Please supply either a uvfits file to interpolate onto (uvfile), a list of uv points (uu & vv), or a gcf_holder cache"



    ###########################
    #   uu and vv retrieval   #
    ###########################

    # if we don't have uu and vv specified, then read them in from the data file
    if (uvfile != 0):
        if "fits" in uvfile:
            data_vis = import_data_uvfits(uvfile)
        elif "ms" in uvfile:
            data_vis = import_data_ms(uvfile)
        else:
            print "not a valid data file to interpolate onto"
        if (verbose==True): print "Read data file to interpolate onto: "+uvfile

    # or read from the cache being fed in
    elif (gcf_holder != 0):
        gcf_holder = gcf_holder
        


    ####################################################
    #   Read the model and apply correction function   #
    ####################################################

    # now that we have either a data_vis or a list of uu,vv points, let's import the model file
    mod_sky_img = import_model(imagefile)
    if (verbose==True): print "Read model file to be interpolated: "+imagefile

    # now apply the correction function
    t0 = time.time()
    if (verbose==True): print "Applying corrfun"

    corr_cache = apply_corrfun(mod_sky_img, mu_RA, mu_DEC, corr_cache=corr_cache, return_cache=return_corr_cache)

    t1 = time.time()
    if (verbose==True): print "corr_fun apply time = " + str(t1-t0)




    ####################################
    #   FFT the image and phase shift  #
    ####################################

    if (verbose==True): print "Starting FFT"
    t0 = time.time()

    mod_fft = transform(mod_sky_img)

    t1 = time.time()
    if (verbose==True): print "fft time = " + str(t1-t0)
    if (verbose==True): print "FFT complete"

    if ((mu_RA != 0) or (mu_DEC != 0)):
        if (verbose==True): print "Starting phase shift"
        t0 = time.time()
        
        phase_shift(mod_fft, mu_RA, mu_DEC)
        
        t1 = time.time()
        if (verbose==True): print "Phase shift time = " + str(t1-t0)
        if (verbose==True): print "Phase shift complete"

    if (verbose==True): print "Starting interpolation"


    ###################################
    #   Do the actual interpolation   #
    ###################################

    t0 = time.time()

    if (uvfile != 0):
        interp, gcf_holder = interpolate_uv(data_vis.uu, data_vis.vv, mod_fft, return_gcf=return_gcf)

    elif (gcf_holder != 0):
        interp, dummy = interpolate_uv(gcf_holder.uu, gcf_holder.vv, mod_fft, gcf_holder=gcf_holder)

    else:
        interp, gcf_holder = interpolate_uv(uu, vv, mod_fft, return_gcf=return_gcf)

    t1 = time.time()
    if (verbose==True): print "interpolation time = " + str(t1-t0)
    if (verbose==True): print "Interpolation complete"


    ########################
    #   Now return stuff   #
    ########################

    # simplest case is just writing to a file:
    if (writefile==True):
        if (verbose==True): print "Writing out to file: "+outfile
        interp_vis = Visibility(interp, data_vis.uu, data_vis.vv, np.ones(interp.shape), data_vis.freqs)
        if "fits" in outfile:
            export_uvfits_from_clone(interp_vis, outfile, uvfile)
        if "ms" in outfile:
            export_ms_from_clone(interp_vis, outfile, uvfile)

        # and we're done!
        return ""

    if (return_gcf == True):
        if (return_corr_cache == True):
            return interp, gcf_holder, corr_cache
        else:
            return interp, gcf_holder
    elif (return_corr_cache == True):
        return interp, corr_cache
    else:
        return interp
