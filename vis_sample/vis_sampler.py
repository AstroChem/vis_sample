import numpy as np
from gridding import *
from classes import *
from transforms import *
from interpolation import interpolate_uv
from file_handling import *
import time

def vis_sample(imagefile, uvfile=0, uu=0, vv=0, writefile=False, outfile="", verbose=True, return_cache=False):

# Error cases #

    if (writefile==True and outfile==""):
        print "Please supply an output imagefile (outfile) if writing to a file"
        return ""

    if (writefile==True and return_cache==True):
        print "Can only return convolution function cache when not writing out to a file"
        return ""




# Either we read in uu and vv coordinates
    if (uvfile == 0):
        if (np.size(uu) == 0) or (np.size(vv) == 0):
            print "Please supply either a uvfits file to interpolate onto (uvfile) or a list of uv points (uu & vv)"
        else:
            mod_sky_img = import_model(imagefile)
            if (verbose==True): print "Read model file to be interpolated: "+imagefile

            if (verbose==True): print "Applying corrfun"
            t0 = time.time()
            apply_corrfun(mod_sky_img, 0.0, 0.0)
            t1 = time.time()
            print "corr_fun apply time = " + str(t1-t0)


            if (verbose==True): print "Starting FFT"
            t0 = time.time()
            mod_fft = transform(mod_sky_img)
            t1 = time.time()
            print "fft time = " + str(t1-t0)

            if (verbose==True): print "FFT complete, starting interpolation"
            t0 = time.time()
            interp = interpolate_uv(uu, vv, mod_fft)
            t1 = time.time()
            print "interpolation time = " + str(t1-t0)
            if (verbose==True): print "Interpolation complete, returning the interpolated visibilities"

            if (writefile==True):
                print "Cannot write file without an input uvfits file (uvfile)"

            return interp




# or we get them from a uvfits file
    else:
        mod_sky_img = import_model(imagefile)
        if (verbose==True): print "Read model file to be interpolated: "+imagefile

        data_vis, data_hd = import_data(uvfile)
        if (verbose==True): print "Read data file to interpolate onto: "+uvfile

        if (verbose==True): print "Applying corrfun"
        t0 = time.time()
        apply_corrfun(mod_sky_img, 0.0, 0.0)
        t1 = time.time()
        print "corr_fun apply time = " + str(t1-t0)

        if (verbose==True): print "Starting FFT"
        t0 = time.time()
        mod_fft = transform(mod_sky_img)
        t1 = time.time()
        print "fft time = " + str(t1-t0)

        if (verbose==True): print "FFT complete, starting interpolation"
        t0 = time.time()
        interp = interpolate_uv(data_vis.uu, data_vis.vv, mod_fft)
        t1 = time.time()
        print "interpolation time = " + str(t1-t0)

        if (verbose==True): print "Interpolation complete"

        if (writefile):
            if (verbose==True): print "Writing out to file: "+outfile
            interp_vis = Visibility(interp, data_vis.uu, data_vis.vv, np.ones(interp.shape), data_vis.lams)
            if "hdf5" in outfile:
                export_hdf5(interp_vis, outfile)
            if "fits" in outfile:
                export_uvfits_from_clone(interp_vis, outfile, uvfile)
        else:
            if (verbose==True): print "Returning interpolated visibilities"
            return interp
    return ""
