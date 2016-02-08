#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="Convert UVHDF5 files into UVFITS files.")
parser.add_argument("HDF5", default="model.hdf5", help="The name of the UVHDF5 file you wish to import.")
parser.add_argument("FITS", help="The original FITS data set, so that we can copy it to stuff in new values.")
parser.add_argument("--out", default="model.vis.fits", help="The output FITS file.")

args = parser.parse_args()

from astropy.io import fits
import h5py
import numpy as np
import shutil

cc = 2.99792458e10 # [cm s^-1]

# Open the old file so we can determine if it had the frequencies increasing or decreasing
f = fits.open(args.FITS)
data = f[0].data
hdr = f[0].header
nfreq = hdr["NAXIS4"]
dnu = hdr["CDELT4"]
ofreqs = hdr["CRVAL4"] + dnu * np.arange(nfreq)  # [Hz] The original frequency values, so we can check that we are inserting in the right order.
f.close()

# Copy the original dataset to something that we can overwrite with the new data set.
shutil.copy(args.FITS, args.out)

# Read the model from the HDF5 file
fid = h5py.File(args.HDF5, "r")
freqs = fid["freqs"][:] # [Hz]
uu = fid["uu"][:,:] # [kilolam]
vv = fid["vv"][:,:] # [kilolam]
real = fid["real"][:,:] # [Jy]
imag = fid["imag"][:,:] # [Jy]
weight = fid["weight"][:,:] #[1/Jy^2]
fid.close()

# Make sure that the original dataset and our UVHDF5 file have the same number of channels and visibilities
assert nfreq == len(freqs), "Number of frequencies mismatched between UVFITS ({}) and UVHDF5 ({}).".format(nfreq, len(freqs))
assert data["DATA"].shape[0] == uu.shape[1], "Number of visibilities mismatched between UVFITS ({}) and UVHDF5 ({}).".format(data["DATA"].shape[0], uu.shape[1])

# Overwrite a copy of the original dataset with these values.
hdulist = fits.open(args.out, mode="update")

# data["DATA"] is originally (nvis, 1, 1, nchan, 1, 3)
# The final dimension (3) is (real, imag, weight)

# Concatenate the different parts of the visibility
D = np.array([real, imag, weight]).T

# Figure out if the original dataset was in reverse order or not
if dnu < 0:
    print("Originally stored frequencies in decreasing order, flipping UVHDF5 to match.")
    D = D[:, ::-1, :]
    freqs = freqs[::-1]
else:
    print("Originally stored with frequencies in increasing order, keeping UVHDF5 order the same.")

assert np.allclose(ofreqs, freqs), "Frequencies between original UVFITS and imported UVHDF5 do not match.\n ofreqs: ({}) \n freqs: ({})".format(ofreqs, freqs)

# Add the zombie dimensions back in to the visibilities
vis = D[:, np.newaxis, np.newaxis, :, np.newaxis, :]

hdulist[0].data["DATA"][:] = vis

# Write the changes to disk
hdulist.flush()
hdulist.close()
