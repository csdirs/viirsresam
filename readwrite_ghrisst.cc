#include "viirsresam.h"

enum {
	MAXDIMS = 5,
};

// Print out error for NetCDF error number n and exit the program.
void
ncfatal(int n, const char *fmt, ...)
{
	va_list args;

	fflush(stdout);
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);

	fprintf(stderr, ": %s\n", nc_strerror(n));
	exit(2);
}

float
ghrsst_readattr(int ncid, int varid, const char *name)
{
	float value;
	nc_type xtype;
	size_t len;
	int n;
	
	n = nc_inq_att(ncid, varid, "add_offset", &xtype, &len);
	if(n != NC_NOERR){
		ncfatal(n, "nc_inq_att failed for add_offset");
	}
	if(xtype != NC_FLOAT || len != 1){
		eprintf("unsupported attribute type/length\n");
	}
	n = nc_get_att_float(ncid, varid, name, &value);
	if(n != NC_NOERR){
		ncfatal(n, "nc_get_att_float failed");
	}
	return value;
}

int
ghrsst_readvar(int ncid, const char *name, Mat &img)
{
	int i, varid, n, ndims, dimids[MAXDIMS], ishape[MAXDIMS], cvt;
	size_t shape[MAXDIMS];
	nc_type nct;
	
	n = nc_inq_varid(ncid, name, &varid);
	if(n != NC_NOERR)
		ncfatal(n, "nc_inq_varid failed for variable %s", name);

	n = nc_inq_var(ncid, varid, NULL, &nct, &ndims, dimids, NULL);
	if(n != NC_NOERR)
		ncfatal(n, "nc_inq_var failed for variable %s", name);
	if(ndims > MAXDIMS)
		eprintf("number of dimensions %d > MAXDIMS=%d\n", ndims, MAXDIMS);
	
	for(i = 0; i < ndims; i++){
		n = nc_inq_dimlen(ncid, dimids[i], &shape[i]);
		if(n != NC_NOERR)
			ncfatal(n, "nc_inq_dimlen failed for dim %d", dimids[i]);
	}
	
	cvt = -1;
	switch(nct){
	default:
		eprintf("unknown netcdf data type");
		break;
	case NC_BYTE:	cvt = CV_8SC1; break;
	case NC_UBYTE:	cvt = CV_8UC1; break;
	case NC_SHORT:	cvt = CV_16SC1; break;
	case NC_USHORT:	cvt = CV_16UC1; break;
	case NC_INT:	cvt = CV_32SC1; break;
	case NC_FLOAT:	cvt = CV_32FC1; break;
	case NC_DOUBLE:	cvt = CV_64FC1; break;
	}
	
	if(ndims == 3 && shape[0] == 1){
		for(i = 1; i < ndims; i++)
			ishape[i-1] = shape[i];
		ndims = 2;
	}else{
		for(i = 0; i < ndims; i++)
			ishape[i] = shape[i];
	}
	
	img.create(ndims, ishape, cvt);
	n = nc_get_var(ncid, varid, img.data);
	if(n != NC_NOERR)
		ncfatal(n, "readvar: nc_get_var failed");
	return varid;
}
