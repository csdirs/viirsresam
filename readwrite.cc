#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// This subroutine reads/writes VIIRS data from/to HDF5 file as unsigned short (2 byte unsigned integer)
// type array, to be scaled with gain and offset values.
//
// Arguments:
// unsigned short **    buffer      IN/OUT    If readwrite == 0, on output contains the pointer to an
//                                            newly allocated array with read data; can be undefined on input;
//                                            but caller is responsible for release of this memory later.
//                                            If readwrite != 0, on input must contain the pointer to
//                                            array with data to be written; unchanged on return.
//
// unsigned long long * dimsizes    OUT       On output contains the sizes of dimensions of 2d array
//                                            dimsizes[0] = swath width (should be 3200, VIIRS)
//                                            dimsizes[1] = height in pixels if granule (depends on size)
//
// float *              gain        OUT       Returns gain parameter for scaling data
// float *              offset      OUT       Returns offset parameter for scaling data
//                                            Physical BT value = gain*buffer[0][i] + offset
//
// char *               filename    IN        Name of HDF5 file from/to which read/write
//
// char *               BTstr       IN        Name of data field in HDF5 file from/to which read/write
//                                            Such as "All_Data/VIIRS-M12-SDR_All/BrightnessTemperature"
//
// int                  readwrite   IN        if readwrite == 0, read data
//                                            if readwrite != 0, write data
//
// Return value:
// Upon sucessful completion, the return value is 0; nonzero return value indicates error.
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int readwrite_viirs(unsigned short **buffer, unsigned long long * dimsizes, float * gain, float * offset,
                    char * filename, char * BTstr, int readwrite)
{

	hid_t   file_id, dataset, dataset_factors, dataspace;
	herr_t  hdferr;
	int     rank_BT, iprint = 0;
	float   gain_offset[2];
	unsigned long long   maxdimsizes[2];

	char BTFstr[256];
	sprintf(BTFstr,"%sFactors", BTstr);

	if(iprint>0) printf("BTstr  = %s\n", BTstr);
	if(iprint>0) printf("BTFstr = %s\n", BTFstr);

	hdferr = H5open();
	if(hdferr!=0) {
		printf("Cannot initialize HDF5 library!\n");
		return -1;
	}

	if(readwrite==0)  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	else              file_id = H5Fopen(filename, H5F_ACC_RDWR,   H5P_DEFAULT);
	if(file_id<0) {
		printf("Cannot open HDF5 file %s!\n", filename);
		return -1;
	}

	dataset_factors = H5Dopen(file_id, BTFstr, H5P_DEFAULT);
	if(dataset_factors<0) {
		printf("Cannot open HDF5 dataset %s!\n", BTFstr);
		return -1;
	}

	hdferr = H5Dread(dataset_factors, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gain_offset);
	if(hdferr<0) {
		printf("Cannot read BT factors!\n");
		return -1;
	}

	*gain   = gain_offset[0];
	*offset = gain_offset[1];
	if(iprint>0) printf("gain = %3.8e  offset = %3.8e\n", *gain, *offset);

	hdferr = H5Dclose(dataset_factors);
	if(hdferr<0) {
		printf("Cannot close HDF5 dataset %s!\n", BTFstr);
		return -1;
	}

	dataset = H5Dopen(file_id, BTstr, H5P_DEFAULT);
	if(dataset<0) {
		printf("Cannot open HDF5 dataset %s!\n", BTstr);
		return -1;
	}

	dataspace = H5Dget_space(dataset);
	if(dataspace<0) {
		printf("Cannot open HDF5 dataspace for dataset %s!\n", BTstr);
		return -1;
	}

	rank_BT = H5Sget_simple_extent_ndims(dataspace);
	if(rank_BT!=2) {
		printf("Unexpected rank of dataspace %i expected 2\n", rank_BT);
		return -1;
	}

	rank_BT = H5Sget_simple_extent_dims(dataspace, dimsizes, maxdimsizes);
	if(rank_BT<0)  {
		printf("Cannot get dataspace dimensions!\n");
		return -1;
	}

	if(iprint>0) printf("dimsizes = %llu  %llu\n", dimsizes[0], dimsizes[1]);


	if(readwrite==0) {
		*buffer = (unsigned short *) malloc(dimsizes[0]*dimsizes[1]*sizeof(unsigned short));
		if(*buffer==NULL) {
			printf("Cannot allocate memory\n");
			return -1;
		}

		hdferr = H5Dread( dataset, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, *buffer);
		if(hdferr<0) {
			printf("Cannot read data to hdf!\n");
			return -1;
		}
	} else {
		hdferr = H5Dwrite(dataset, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, *buffer);
		if(hdferr<0) {
			printf("Cannot write data to hdf!\n");
			return -1;
		}
	}

	// close everything
	hdferr = H5Dclose(dataset);
	if(hdferr<0) {
		printf("Cannot close HDF5 dataset %s!\n", BTstr);
		return -1;
	}

	hdferr = H5Fclose(file_id);
	if(hdferr<0) {
		printf("Cannot close HDF5 file %s!\n", filename);
		return -1;
	}

	hdferr = H5close();
	if(hdferr<0) {
		printf("Cannot close HDF5 library!\n");
		return -1;
	}

	return 0;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// This subroutine reads/writes VIIRS data from/to HDF5 file as float (4 byte floating point)
// type array. Used for VIIRS band M13.
//
// Arguments:
// float **             buffer      IN/OUT    If readwrite = 0, on output contains the pointer to an
//                                            newly allocated array with read data; can be undefined on input.
//                                            If readwrite != 0, on input must contain the pointer to
//                                            array with data to be written; unchanged on return.
//
// unsigned long long * dimsizes    OUT       On output contains the sizes of dimensions of 2d array
//                                            dimsizes[0] = swath width (should be 3200, VIIRS)
//                                            dimsizes[1] = height in pixels if granule (depends on size)
//
// char *               filename    IN        Name of HDF5 file from/to which read/write
//
// char *               BTstr       IN        Name of data field in HDF5 file from/to which read/write
//                                            Such as "All_Data/VIIRS-M12-SDR_All/BrightnessTemperature"
//
// int                  readwrite   IN        if readwrite = 0, read data
//                                            if readwrite != 0, write data
//
// Return value:
// Upon sucessful completion, the return value is 0; nonzero return value indicates error.
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int readwrite_viirs_float(float **buffer, unsigned long long * dimsizes, char * filename, char * BTstr, int readwrite)
{

	hid_t   file_id, dataset, dataspace;
	herr_t  hdferr;
	int     rank_BT, iprint = 0;
	unsigned long long   maxdimsizes[2];

	if(iprint>0) printf("BTstr  = %s\n", BTstr);

	hdferr = H5open();
	if(hdferr!=0) {
		printf("Cannot initialize HDF5 library!\n");
		return -1;
	}

	if(readwrite==0)  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	else              file_id = H5Fopen(filename, H5F_ACC_RDWR,   H5P_DEFAULT);
	if(file_id<0) {
		printf("Cannot open HDF5 file %s!\n", filename);
		return -1;
	}

	dataset = H5Dopen(file_id, BTstr, H5P_DEFAULT);
	if(dataset<0) {
		printf("Cannot open HDF5 dataset %s!\n", BTstr);
		return -1;
	}

	dataspace = H5Dget_space(dataset);
	if(dataspace<0) {
		printf("Cannot open HDF5 dataspace for dataset %s!\n", BTstr);
		return -1;
	}

	rank_BT = H5Sget_simple_extent_ndims(dataspace);
	if(rank_BT!=2) {
		printf("Unexpected rank of dataspace %i expected 2\n", rank_BT);
		return -1;
	}

	rank_BT = H5Sget_simple_extent_dims(dataspace, dimsizes, maxdimsizes);
	if(rank_BT<0)  {
		printf("Cannot get dataspace dimensions!\n");
		return -1;
	}

	if(iprint>0) printf("dimsizes = %llu  %llu\n", dimsizes[0], dimsizes[1]);



	if(readwrite==0) {

		*buffer = (float *) malloc(dimsizes[0]*dimsizes[1]*sizeof(float));
		if(*buffer==NULL) {
			printf("Cannot allocate memory\n");
			return -1;
		}

		hdferr = H5Dread( dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, *buffer);
		if(hdferr<0) {
			printf("Cannot read data from hdf!\n");
			return -1;
		}
	} else {

		hdferr = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, *buffer);
		if(hdferr<0) {
			printf("Cannot write data to hdf!\n");
			return -1;
		}

	}

	// close everything
	hdferr = H5Dclose(dataset);
	if(hdferr<0) {
		printf("Cannot close HDF5 dataset %s!\n", BTstr);
		return -1;
	}

	hdferr = H5Fclose(file_id);
	if(hdferr<0) {
		printf("Cannot close HDF5 file %s!\n", filename);
		return -1;
	}

	hdferr = H5close();
	if(hdferr<0) {
		printf("Cannot close HDF5 library!\n");
		return -1;
	}

	return 0;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This subroutine writes an attribute for a particular field in a HDF5 file.
//
// Arguments:
//
// char *    filename      IN      Name of HDF5 file to which write an attribute
//
// char *    attrFieldStr  IN      Name of data field in HDF5 file to which write an attribute
//                                 Such as "All_Data/VIIRS-M12-SDR_All/BrightnessTemperature"
//
// char *    attrNameStr   IN      Name of the attribute to write, such as "Resampling"
//
// float     destrval      IN      Value of the attribute to write
//
// Return value:
// Upon sucessful completion, the return value is nonnegative;
//                negative return value indicates error;
//                positive return value indicates that the attribute was already set.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int write_viirs_attribute(char * filename, char * attrFieldStr, char * attrNameStr, float destrval)
{

	hid_t   file_id, dataset;
	herr_t  hdferr;
	int     retval = 0;

	hdferr = H5open();
	if(hdferr!=0) {
		printf("Cannot initialize HDF5 library!\n");
		return -1;
	}

	file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if(file_id<0) {
		printf("Cannot open HDF5 file %s!\n", filename);
		return -1;
	}

	dataset = H5Dopen(file_id, attrFieldStr, H5P_DEFAULT);
	if(dataset<0) {
		printf("Cannot open HDF5 dataset %s!\n", attrFieldStr);
		return -1;
	}

	if(H5Aexists(dataset, attrNameStr)>0) {
		// attribute already exists
		retval = 1;
	} else {
		// attribute does not exist, create it
		hid_t type_id = H5Tcopy(H5T_NATIVE_FLOAT);
		if(type_id<0)  {
			printf("Cannot create a new datatype!\n");
			return -1;
		}

		hid_t space_id = H5Screate(H5S_SCALAR);
		if(space_id<0)  {
			printf("Cannot create a new dataspace!\n");
			return -1;
		}

		hid_t attr_id = H5Acreate( dataset, attrNameStr, type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
		if(attr_id<0) {
			printf("Cannot create a HDF5 attribute\n");
			return -1;
		}

		float destrval = 1.0;
		float *pdestrval = &destrval;
		hdferr = H5Awrite( attr_id, H5T_NATIVE_FLOAT, (const void *) pdestrval );
		if(hdferr==-1) {
			printf("Cannot write attribute!\n");
			return -1;
		}
	}

	// close everything
	hdferr = H5Dclose(dataset);
	if(hdferr<0) {
		printf("Cannot close HDF5 dataset %s!\n", attrFieldStr);
		return -1;
	}

	hdferr = H5Fclose(file_id);
	if(hdferr<0) {
		printf("Cannot close HDF5 file %s!\n", filename);
		return -1;
	}

	hdferr = H5close();
	if(hdferr<0) {
		printf("Cannot close HDF5 library!\n");
		return -1;
	}

	return retval;
};
