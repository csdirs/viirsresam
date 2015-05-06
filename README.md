## VIIRS Resampling

This program resamples VIIRS data.

Dependencies:

* C++ toolchain
* OpenCV
* HDF5 library

Run `make` to build the program named `viirsresam`. Running the program
on a granule will modify the data in-place and add an attribute indicating
it was resampled.
