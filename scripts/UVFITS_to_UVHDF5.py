#!/usr/bin/env python

# This is designed as a reference implementation of the UVHDF5 specification. It shows a conversion of SMA data in the UVFITS format to the UVHDF5 format.

import argparse

parser = argparse.ArgumentParser(description="Convert SMA FITS files into UVHDF5 file format.")
parser.add_argument("FITS", help="The input UVFITS file.")
parser.add_argument("--out", default="data.hdf5", help="The output UVHDF5 file.")
args = parser.parse_args()

from astropy.io import fits
import h5py
import numpy as np

cc = 2.99792458e10 # [cm s^-1]

# Reading SMA dataset

f = fits.open(args.FITS)

data = f[0].data
hdr = f[0].header
nfreq = hdr["NAXIS4"]
dnu = hdr["CDELT4"]

freqs = hdr["CRVAL4"] + dnu * np.arange(nfreq)  # Hz

# Convert each uu, vv cordinate from light nanoseconds to kilolambda,
# depending on which channel frequency we are using
uu = 1e-3 * (freqs * np.tile(data["UU"], (nfreq, 1)).T).T
vv = 1e-3 * (freqs * np.tile(data["VV"], (nfreq, 1)).T).T

# uu, vv are now (nfreq, nvis) shape arrays

shape = uu.shape
nvis = uu.shape[1]

# print("Original shape of DATA", data["DATA"].shape)

# Remove all of the "zombie" 1D columns
vis = np.squeeze(data["DATA"])  

# Now, vis is stored as an (npoints, nfreqs, 3) array, where last dimension is
# (real, imag, weight)

# Read and convert all of these to (nfreq, nvis) arrays
real = vis[:, :, 0].T
imag = vis[:, :, 1].T
weight = vis[:, :, 2].T

# Now, stuff each of these into an HDF5 file.
    fid = h5py.File(args.out, "w")

# Add in observational attributes
for key in ["OBJECT", "TELESCOP", "ORIGIN"]:
    try:
        val = hdr[key]
        fid.attrs[key] = val
    except KeyError:
        continue

# Add in format specification version
fid.attrs["FMT_Version"] = "v0.1"

# Are the frequencies stored in increasing or decreasing order in UVFITS?
# UVHDF5 always stores frequencies in increasing order
if dnu > 0:
    fid.create_dataset("freqs", (nfreq,), dtype="float64")[:] = freqs # [Hz]

    fid.create_dataset("uu", shape, dtype="float64")[:,:] = uu # [kilolambda]
    fid.create_dataset("vv", shape, dtype="float64")[:,:] = vv # [kilolambda]

    fid.create_dataset("real", shape, dtype="float64")[:,:] = real # [Jy]
    fid.create_dataset("imag", shape, dtype="float64")[:,:] = imag # [Jy]

    fid.create_dataset("weight", shape, dtype="float64")[:,:] = weight #[1/Jy^2]

else:
    print("UVFITS stored frequencies in decreasing order, flipping to positive for UVHDF5")
    fid.create_dataset("freqs", (nfreq,), dtype="float64")[:] = freqs[::-1] # [Hz]

    fid.create_dataset("uu", shape, dtype="float64")[:,:] = uu[::-1] # [kilolambda]
    fid.create_dataset("vv", shape, dtype="float64")[:,:] = vv[::-1] # [kilolambda]

    fid.create_dataset("real", shape, dtype="float64")[:,:] = real[::-1] # [Jy]
    fid.create_dataset("imag", shape, dtype="float64")[:,:] = imag[::-1] # [Jy]

    fid.create_dataset("weight", shape, dtype="float64")[:,:] = weight[::-1] #[1/Jy^2]

fid.close()
