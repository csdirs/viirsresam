#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include "viirsresam.h"

#define TRUE	1
#define FALSE	0

// Write data in HDF5 file named filename with layer named varname.
// The layout will be created if it doesn't exist already.
//
void
create_viirs(Mat data, const char *filename, const char *varname)
{
	hid_t dataset, dataspace, dtype;
	
	CV_Assert(data.dims == 2);
	switch(data.type()){
	default:
		eprintf("unsupported Mat type %d\n", data.type());
	case CV_16UC1:
		dtype = H5T_NATIVE_USHORT;
		break;
	case CV_32FC1:
		dtype = H5T_NATIVE_FLOAT;
		break;
	}
	
	herr_t  hdferr = H5open();
	if(hdferr < 0){
		eprintf("cannot initialize HDF5 library:");
	}
	hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if(file_id < 0){
		eprintf("cannot open HDF5 file %s", filename);
	}
	
	if(H5Lexists(file_id, varname, H5P_DEFAULT) == TRUE){
		// dataset exists
		hsize_t dims[2], maxdims[2];
		
		dataset = H5Dopen(file_id, varname, H5P_DEFAULT);
		if(dataset < 0){
			eprintf("cannot open HDF5 dataset %s", varname);
		}
		dataspace = H5Dget_space(dataset);
		if(dataspace<0){
			eprintf("cannot open HDF5 dataspace for dataset %s", varname);
		}
		int rank = H5Sget_simple_extent_ndims(dataspace);
		if(rank != 2){
			eprintf("unexpected rank %d of dataspace; expected 2", rank);
		}
		rank = H5Sget_simple_extent_dims(dataspace, dims, maxdims);
		if(rank != 2 || (int)dims[0] != data.rows || (int)dims[1] != data.cols){
			eprintf("HDF5 dataspace dimensions is %dx%d; expected %dx%d",
				dims[0], dims[1], data.rows, data.cols);
		}
	}else{
		// dataset does not exist, so create it
		hsize_t dims[2];
		dims[0] = data.rows;
		dims[1] = data.cols;
		
		dataspace = H5Screate_simple(2, dims, NULL);
		if(dataspace < 0){
			eprintf("cannot create HDF5 dataspace for dataset %s", varname);
		}
		dataset = H5Dcreate(file_id, varname, dtype, dataspace, 
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(dataset < 0){
			eprintf("cannot create HDF5 dataset %s", varname);
		}
	}

	hdferr = H5Dwrite(dataset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data);
	if(hdferr < 0){
		eprintf("Cannot write data to hdf");
	}

	// close everything
	hdferr = H5Dclose(dataset);
	if(hdferr < 0){
		eprintf("cannot close HDF5 dataset %s", varname);
	}
	hdferr = H5Sclose(dataspace);
	if(hdferr < 0){
		eprintf("cannot close HDF5 dataset %s", varname);
	}
	hdferr = H5Fclose(file_id);
	if(hdferr < 0){
		eprintf("cannot close HDF5 file %s", filename);
	}
	hdferr = H5close();
	if(hdferr < 0){
		eprintf("cannot close HDF5 library");
	}
}
