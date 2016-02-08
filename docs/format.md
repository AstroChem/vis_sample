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


