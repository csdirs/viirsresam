// VIIRS resampling

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "viirsresam.h"

enum {
	UNKNOWN,
	L2P_GHRSST,
	GMODO,
	GMTCO,
};

#define _LATNAME	"All_Data/VIIRS-MOD-GEO_All/Latitude"
#define _LONNAME	"All_Data/VIIRS-MOD-GEO_All/Longitude"
#define _TCLATNAME	"All_Data/VIIRS-MOD-GEO-TC_All/Latitude"
#define _TCLONNAME	"All_Data/VIIRS-MOD-GEO-TC_All/Longitude"
#define LATNAME _LATNAME
#define LONNAME _LONNAME
#define GEO_RESAM_ATTR_NAME	"Resampling"

#define GETARG(x)	do{\
		(x) = *argv++;\
		argc--;\
	}while(0);

char *progname;

static void
usage()
{
	printf("usage: %s [-s] geofile viirs_h5_file\n", progname);
	printf("       %s geofile\n", progname);
	printf("\n");
	printf("	-s	reorder output by latitude\n");
	printf("\n");
	printf("If viirs_h5_file is given, brightness temperature is resampled\n");
	printf("and saved in the input file, in addition to a \"Resampling\"\n");
	printf("attribute indicating the data is already resampled.\n");
	printf("\n");
	printf("If viirs_h5_file is not given, the latitude and longitude in geofile\n");
	printf("is reordered to the same order resulting from the -s flag. Also,\n");
	printf("a \"Resampling\" attribute is added indicating latitude/longitude\n");
	printf("is already reordered.\n");
	exit(2);
}

static char*
filebasename(char *path)
{
	char *p = strrchr(path, '/');
	if(p == NULL){
		return path;
	}
	return p+1;
}

static int
getfiletype(char *path)
{
	char *p = filebasename(path);
	if(strlen(p) > 20 && strncmp(&p[20], "L2P_GHRSST", strlen("L2P_GHRSST")) == 0){
		return L2P_GHRSST;
	}
	if(strlen(p) > 20 && strncmp(&p[20], "L2P_GHRSST", strlen("L2P_GHRSST")) == 0){
		return L2P_GHRSST;
	}
	if(strncmp(p, "GMODO_npp_", strlen("GMODO_npp_")) == 0){
		return GMODO;
	}
	if(strncmp(p, "GMTCO_npp_", strlen("GMTCO_npp_")) == 0){
		return GMTCO;
	}
	return UNKNOWN;
}

static void
sortlatlon(const char *geofile)
{
	Mat sind;
	int status;
	unsigned long long dims1[32];
	float *latbuf = NULL;
	float *lonbuf = NULL;
	
	// read latitude & longitude
	status = readwrite_viirs_float(&latbuf, dims1, geofile, LATNAME, 0);
	if(status != 0){
		fprintf(stderr, "Cannot read VIIRS (lat) geolocation data!\n");
		exit(2);
	}
	status = readwrite_viirs_float(&lonbuf, dims1, geofile, LONNAME, 0);
	if(status != 0){
		fprintf(stderr, "Cannot read VIIRS (lon) geolocation data!\n");
		exit(2);
	}
	int sy = dims1[0];	// height, along the track
	int sx = dims1[1];	// width, across track, along scan line
	Mat lat(sy, sx, CV_32FC1, latbuf);
	Mat lon(sy, sx, CV_32FC1, lonbuf);

	// sort latitude & longitude
	getsortingind(sind, sy);
	Mat slat = resample_sort(sind, lat);
	Mat slon = resample_sort(sind, lon);
	CHECKMAT(slat, CV_32FC1);
	CHECKMAT(slon, CV_32FC1);
	
	// write sorted latitude & longitude
	status = readwrite_viirs_float((float**)&slat.data, dims1, geofile, LATNAME, 1);
	if(status != 0){
		fprintf(stderr, "Cannot read VIIRS (lat) geolocation data!\n");
		exit(2);
	}
	status = readwrite_viirs_float((float**)&slon.data, dims1, geofile, LONNAME, 1);
	if(status != 0){
		fprintf(stderr, "Cannot read VIIRS (lon) geolocation data!\n");
		exit(2);
	}

	// write resampling attribute for latitude & longitude
	int estat = 0;
	status = write_viirs_attribute(geofile, LATNAME, GEO_RESAM_ATTR_NAME, 1.0);
	if(status < 0){
		printf("ERROR: Cannot write VIIRS attribute!\n");
		estat = 2;
	}
	if(status > 0){
		printf("WARNING! Data was already resampled\n");
	}
	status = write_viirs_attribute(geofile, LONNAME, GEO_RESAM_ATTR_NAME, 1.0);
	if(status < 0){
		printf("ERROR: Cannot write VIIRS attribute!\n");
		estat = 2;
	}
	if(status > 0){
		printf("WARNING! Data was already resampled\n");
	}

	free(latbuf);
	free(lonbuf);
	exit(estat);
}

static void
run_ghrsst(char *ncfile, bool sortoutput)
{
	int ncid, n;
	Mat _sst, lat, lon;
	
	n = nc_open(ncfile, NC_WRITE, &ncid);
	if(n != NC_NOERR)
		ncfatal(n, "nc_open failed for %s", ncfile);
	
	int varid = ghrsst_readvar(ncid, "sea_surface_temperature", _sst);
	float offset = ghrsst_readattr(ncid, varid, "add_offset");
	float scale = ghrsst_readattr(ncid, varid, "scale_factor");
	printf("scale = %f, offset = %f\n", scale, offset);
	
	CHECKMAT(_sst, CV_16SC1);
	short *sst = (short*)_sst.data;
	Mat _sstf = Mat::zeros(_sst.size(), CV_32FC1);
	float *sstf = (float*)_sstf.data;
	for(int i = 0; i < (int)_sst.total(); i++){
		sstf[i] = sst[i]*scale + offset;
	}
	
	ghrsst_readvar(ncid, "lat", lat);
	ghrsst_readvar(ncid, "lon", lon);
	resample_viirs_mat(_sstf, lat, lon, -999, sortoutput);
}

static void
run_tcgeo(char *gmodofile, char *gmtcofile, bool sortoutput)
{
	int status;
	unsigned long long dims[32];
	float *buflat = NULL;
	float *buflon = NULL;
	float *buftclat = NULL;
	float *buftclon = NULL;
	
	status = readwrite_viirs_float(&buflat, dims, gmodofile, _LATNAME, 0);
	if(status != 0){
		eprintf("Cannot read VIIRS (lat) geolocation data!");
	}
	status = readwrite_viirs_float(&buflon, dims, gmodofile, _LONNAME, 0);
	if(status != 0){
		eprintf("Cannot read VIIRS (lon) geolocation data!\n");
	}
	status = readwrite_viirs_float(&buftclat, dims, gmtcofile, _TCLATNAME, 0);
	if(status != 0){
		eprintf("Cannot read VIIRS (lat) terrain-corrected geolocation data!");
	}
	status = readwrite_viirs_float(&buftclon, dims, gmtcofile, _TCLONNAME, 0);
	if(status != 0){
		eprintf("Cannot read VIIRS (lon) terrain-corrected geolocation data!\n");
	}
	printf("done reading data\n");
	Mat lat(dims[0], dims[1], CV_32FC1, buflat);
	Mat lon(dims[0], dims[1], CV_32FC1, buflon);
	Mat tclat(dims[0], dims[1], CV_32FC1, buftclat);
	Mat tclon(dims[0], dims[1], CV_32FC1, buftclon);
	
	Mat latdiff = tclat - lat;
	Mat londiff = tclon - lon;
	
	printf("resampling lat\n");
	if(DEBUG)dumpmat("tclat.bin", tclat);
	resample_viirs_mat(latdiff, lat, lon, -999, sortoutput);
	Mat tclatp = lat + latdiff;
	if(DEBUG)dumpmat("tclatp.bin", tclatp);
	
	printf("resampling lon\n");
	if(DEBUG)dumpmat("tclon.bin", tclon);
	resample_viirs_mat(londiff, lat, lon, -999, sortoutput);
	Mat tclonp = lon + londiff;
	if(DEBUG)dumpmat("tclonp.bin", tclonp);

	if(DEBUG)exit(3);
	
	free(buflat);
	free(buflon);
	free(buftclat);
	free(buftclon);
}

static int
run_band(char *h5file, char *geofile, bool sortoutput)
{
	unsigned short *buffer1  = NULL;
	float          *bufferf1 = NULL;
	float          *bufferf2 = NULL;
	float          *bufferf3 = NULL;
	unsigned long long dims1[32];
	int status, i, j, is;
	float scale1, offset1;
	int ix, sx, sy;
	double scale, offset;
	float ** img_in, **lat, **lon;

	// extract the name of band from the file
	is = 0;
	for(i=(strlen(h5file)-6); i>=0; i--) {
		// look for sequence "SVM", then following two chars give band number
		if( (h5file[i]=='S') && (h5file[i+1]=='V') && (h5file[i+2]=='M') ) {
			is = (h5file[i+3]-'0')*10 + (h5file[i+4]-'0');
			break;
		}
	}
	printf("Band = %i\n", is);

	// check that the band number extracted from file name is within the valid range
	if((is<1)||(is>16)) {
		printf("ERROR: Invalid band\n");
		return -7;
	}

	// prepare all data field strings for later use
	char attrfieldstr[128], attrnamestr[128], btstr[128];

	// resampling attribute field
	sprintf(attrfieldstr,"Data_Products/VIIRS-M%i-SDR/VIIRS-M%i-SDR_Aggr", is, is);
	printf("Resampling atribute location = %s\n", attrfieldstr);

	// generate names of main data fields to be resampled
	// and the names of the corresponding resampling attributes
	if(is<12) {
		// for M11 and below, resample Reflectance
		sprintf(btstr, "All_Data/VIIRS-M%i-SDR_All/Reflectance", is);
		sprintf(attrnamestr, "ResamplingReflectance");
	} else {

		// for M12 and above, resample Brightness Temperature
		sprintf(btstr, "All_Data/VIIRS-M%i-SDR_All/BrightnessTemperature", is);
		sprintf(attrnamestr, "ResamplingBrightnessTemperature");
	}
	printf("Resampling atribute name = %s\n", attrnamestr);
	printf("Data location = %s\n", btstr);          // name of main data field to be resampled

	// read data
	if(is!=13) {
		status = readwrite_viirs(&buffer1, dims1, &scale1, &offset1, h5file, btstr, 0);
	} else {
		status = readwrite_viirs_float( &bufferf1, dims1, h5file, btstr, 0);
	}
	if(status!=0) {
		printf("ERROR: Cannot read VIIRS data!\n");
		return 10*status;
	}

	// read geolocation data
	status = readwrite_viirs_float( &bufferf2, dims1, geofile, LATNAME, 0);
	if(status!=0) {
		fprintf(stderr, "Cannot read VIIRS (lat) geolocation data!\n");
		exit(1);
	}
	status = readwrite_viirs_float( &bufferf3, dims1, geofile, LONNAME, 0);
	if(status!=0) {
		fprintf(stderr, "Cannot read VIIRS (lon) geolocation data!\n");
		exit(1);
	}

	// extract scale, offset and dimensions info
	sy = dims1[0]; // height, along the track
	sx = dims1[1]; // width, across track, along scan line
	printf("nx = %i ny = %i\n", sx, sy);
	scale = offset = 0;	// quiet compiler warning
	if(is!=13) {
		scale  = ((double) scale1);
		offset = ((double) offset1);
		printf("scale = %f offset = %f\n", scale, offset);
	}

	// allocate temporary data arrays
	img_in = allocate_2d_f(sy, sx);
	if( (img_in==NULL)) {
		printf("ERROR: Cannot allocate memory\n");
		return -1;
	}

	// if needed, apply scale and offset to get physical data
	if(is!=13) {
		// apply scale and offset to get real physical units
		for(ix=0; ix<sx*sy; ix++) {
			img_in[0][ix] = scale*buffer1[ix] + offset;
		}
	} else {
		// no scaling for band 13
		for(ix=0; ix<sx*sy; ix++) {
			img_in[0][ix] = bufferf1[ix];
		}
	}

	// allocate geolocation arrays
	lat      = allocate_2d_f(sy, sx);
	lon      = allocate_2d_f(sy, sx);
	if( (lat==NULL) || (lon==NULL) ) {
		printf("ERROR: Cannot allocate memory\n");
		return -1;
	}

	for(ix=0; ix<sx*sy; ix++) {
		lat[0][ix] = bufferf2[ix];
	}
	for(ix=0; ix<sx*sy; ix++) {
		lon[0][ix] = bufferf3[ix];
	}

	free(bufferf2);
	free(bufferf3);

	// resampling of image on sorted lon, lat grid
	if(is != 13){
		resample_viirs(img_in, lat, lon, sx, sy, scale*DELETION_ZONE_INT+offset, sortoutput);
	}else{
		resample_viirs(img_in, lat, lon, sx, sy, DELETION_ZONE_FLOAT, sortoutput);
	}

	// Scale resampled data back to integers if band != M13
	if(is!=13) {
		for(ix=0; ix<sx*sy; ix++) {

			// scale resampled data back to integer value
			j = (int) round((img_in[0][ix] - offset)/scale);

			// check if integer is in the valid range
			if(j<0) {
				printf("Output data out of range at ( %5i %5i ): %i\n", ix%sx, ix/sx, j);
				j = 0;
			}
			if(j>65535) {
				printf("Output data out of range at ( %5i %5i ): %i\n", ix%sx, ix/sx, j);
				j = 65535;
			}

			buffer1[ix] = (unsigned short) j;

		} // for ix = all data in 2d array

		// write resampled data back to file as short int
		status = readwrite_viirs(&buffer1, dims1, &scale1, &offset1, h5file, btstr, 1);
		free(buffer1);
	} else {
		// no conversion for band M13
		for(ix=0; ix<sx*sy; ix++) {
			bufferf1[ix] = img_in[0][ix];
		} // for ix = all data in 2d array

		// write resampled band M13 data back to file as float
		status = readwrite_viirs_float(&bufferf1, dims1, h5file, btstr, 1);
		free(bufferf1);
	}

	if(status!=0) {
		printf("ERROR: Cannot write VIIRS data!\n");
		return 10*status;
	}

	// write a resampled attribute
	status = write_viirs_attribute(h5file, attrfieldstr, attrnamestr, 1.0);
	if(status < 0){
		printf("ERROR: Cannot write VIIRS attribute!\n");
		return 40*status;
	}
	if(status > 0){
		printf("WARNING! Data was already resampled\n");
	}

	free(img_in[0]);
	free(img_in);
	free(lat[0]);
	free(lat);
	free(lon[0]);
	free(lon);

	return 0;
}

int
main(int argc, char** argv)
{
	char *flag;

	// parse arguments
	GETARG(progname);
	bool sortoutput = false;
	while(argc > 0 && strlen(argv[0]) == 2 && argv[0][0] == '-') {
		GETARG(flag);

		switch(flag[1]) {
		default:
			usage();
			break;
		case '-':
			goto argdone;
		case 's':
			sortoutput = true;
			break;
		}
	}
argdone:
	if(argc == 1 && getfiletype(argv[0]) == L2P_GHRSST){
		printf("resampling GHRSST file...\n");
		run_ghrsst(argv[0], sortoutput);
		exit(0);
	}
	if(argc == 1){
		sortlatlon(argv[0]);
		exit(0);
	}
	if(argc == 2 && getfiletype(argv[0]) == GMODO && getfiletype(argv[1]) == GMTCO){
		run_tcgeo(argv[0], argv[1], sortoutput);
		exit(0);
	}
	if(argc != 2)
		usage();
	char *geofile = argv[0];
	char *h5file = argv[1];

	// echo command line
	printf("viirsresam %s%s %s\n",
		sortoutput ? "-s " : "",
		geofile, h5file);
	printf("Corresponding geofile = %s\n", geofile);

	return run_band(h5file, geofile, sortoutput);
}
