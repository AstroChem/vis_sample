# Format

This document describes the file format used to store the visibilities. It will be updated as this package evolves. For now, it's called the `UVHDF5` format.

Some limitations of this format:

* Only unflagged data will be stored in the HDF5.
* No autocorrelations will be stored in the HDF5.

Package version: v 0.1

Stored in an HDF5 format.

## Header

    Target Name: 
    Creation Date: 
    Reduction Software (name and version): 
    Format Spec Version (vis_sample version):

If it so happens that the arrays are all the same size, there will be an attribute in the header that specifies this. 

    nvis_const = True

## Visibilities

Each channel is stored as a folder within the HDF5, with the name chanX, where X is a sequential index starting at 0. Channels should be stored in strictly increasing frequency, where chan0 is the lowest frequency.

    chan0
    chan1
    chan2
    chan3
    ...
    chan10
    ...
    chan100

Although this extra complexity may seem unnecessary for simple datasets, this is because we want the capacity to easily store datasets that contain channelized flags. If we were to store everything in equal-sized 2D arrays, then we get into the ugly business of needing to propagate flags that the user needs to interpret.

There are required attributes to each folder:

    nu: [Hz]
    dnu: [Hz]

Within the channel folder, we store the `uu` and `vv` points as 1D, length `nvis` arrays. The real and imaginary components of the visibilities are stored in `real` and `imag`, and the weights are stored as `weights`. All of these components should be 1D arrays of length `nvis`, where `nvis` is the number of visibilities in that specific channel, and need not be equal across channels.

If storing spectral line data, you must at least store the visibilities in chan0.

