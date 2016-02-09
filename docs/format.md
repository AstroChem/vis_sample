# UVHDF5 Format

Format specification version: v 0.1

This document describes the file format used to store the visibilities. It will be updated as this package evolves. For now, it's called the `UVHDF5` format.

Some limitations of this format:

* Only unflagged data will be stored in the HDF5.
* No autocorrelations will be stored in the HDF5.

# Rationale of using HDF5 vs. NPZ and other formats

While NPZ is readily readable by many of our scripts already, the HDF5 file format presents many attractive features. 

* HDF5 is a universally supported standard that can be read by programs written in other languages: Python, IDL, MATLAB, C, C++, Fortran, Julia, etc. NPZ is a python-specific format.
* HDF5 supports the use of metadata attributes, much like a FITS header. This means that extra information (e.g., name of the target, original owner of the dataset, ALMA program of observation, reduction version) can be stored inside of the visibility file itself so that we don't get confused about provenance later on. NPZ does not support metadata, so this information can be easily lost if a file is sitting on disk for a few months or is transferred between members of the group.
* Programs like `hdf-java` allow easy visual inspection of the dataset contents (via GUI) to check that everything is in order. Checking a NPZ file field by field can be very tedious through the Python interpreter.

For now, this format has been designed around the most common type of use case: that there are the **same number of visibilities for each channel**. For more complex types of data, we will need to collaborate on a more advanced format.

## Header

The following fields are stored as attributes on the `/` folder of the HDF5 file.

    Target Name: 
    Creation Date: 
    Reduction Software (name and version): 
    Format Spec Version (vis_sample version):

## Frequencies

The frequencies are stored in Hz as a 1D array of length `nchan`. 

    freqs [Hz]

## Visibilities

The visibilities are stored as `(nchan, nvis)` datasets on the HDF5 file.

    uu [kilolambda]
    vv [kilolambda]
    real [Jy]
    imag [Jy]
    weight [1/Jy^2]

# Reading from UVHDF5 

Here is an example snippet of python code that would be useful to read this new file format.

    import h5py
    import numpy as np

    filename = "data.hdf5"

    fid = h5py.File(filename, "r")
    freqs = fid["freqs"][:] # [Hz]
    uu = fid["uu"][:,:] # [kilolam]
    vv = fid["vv"][:,:] # [kilolam]
    real = fid["real"][:,:] # [Jy]
    imag = fid["imag"][:,:] # [Jy]
    weight = fid["weight"][:,:] #[1/Jy^2]
    fid.close()

    # Do fancy stuff with visibilities here


# Writing to UVHDF5

Here is an example of how to write your dataset to this file format in python.

    # Assume that you have your frequencies stored in a 1D array of length `nfreq`
    # and that you have your visibilities stored in 2D arrays of size `(nfreq, nvis)`.

    shape = (nfreq, nvis)

    import h5py

    filename = "model.hdf5"

    fid = h5py.File(filename, "w")
    fid.create_dataset("freqs", (nfreq,), dtype="float64")[:] = freqs # [Hz]

    fid.create_dataset("uu", shape, dtype="float64")[:,:] = uu # [kilolambda]
    fid.create_dataset("vv", shape, dtype="float64")[:,:] = vv # [kilolambda]

    fid.create_dataset("real", shape, dtype="float64")[:,:] = real # [Jy]
    fid.create_dataset("imag", shape, dtype="float64")[:,:] = imag # [Jy]

    fid.create_dataset("weight", shape, dtype="float64")[:,:] = weight #[1/Jy^2]
    fid.close()
