# UVHDF5 Format

This document describes the file format used to store the visibilities. It will be updated as this package evolves. For now, it's called the `UVHDF5` format.

Some limitations of this format:

* Only unflagged data will be stored in the HDF5.
* No autocorrelations will be stored in the HDF5.

Package version: v 0.1

Stored in an HDF5 format.

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
