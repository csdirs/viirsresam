#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <math.h>
#include <opencv2/opencv.hpp>
#include <netcdf.h>

using namespace cv;

typedef unsigned long long uvlong;
typedef unsigned short ushort;

#define nelem(x)	(sizeof(x)/sizeof((x)[0]))
#define CHECKMAT(M, T)	CV_Assert((M).type() == (T) && (M).isContinuous())
#define SQ(x)	((x)*(x))
#define SIGN(A)   ((A) > 0 ? 1 : ((A) < 0 ? -1 : 0 ))
#define RADIANCE(x)	((x)*M_PI/180.0)
#define DEGREE(x)	((x)*180.0/M_PI)

enum {
	VIIRS_WIDTH = 3200,
	NDETECTORS = 16,
	MAX_TEMP = 350,	// in Kelvin
	MIN_TEMP = 0,	// in Kelvin
	INVALID_TEMP = -999,
	DEBUG = false,
	
	DELETION_ZONE_INT = 65533,
	DELETION_ZONE_FLOAT = -999,
};


// allocate_2d.cc
float ** allocate_2d_f(int n1, int n2);
int ** allocate_2d_i(int n1, int n2);

// readwrite.cc
int readwrite_viirs(unsigned short **buffer, unsigned long long * dimsizes, float * gain, float * offset, 
                    char * filename, char * BTstr, int readwrite);
int readwrite_viirs_float(float **buffer, unsigned long long * dimsizes, const char * filename, const char * BTstr, int readwrite);
int write_viirs_attribute(const char *filename, const char *attrFieldStr, const char *attrNameStr, float destrval);

// readwrite_ghrisst.cc
void ncfatal(int n, const char *fmt, ...);
float ghrsst_readattr(int ncid, int varid, const char *name);
int ghrsst_readvar(int ncid, const char *name, Mat &img);

// resample.cc
void resample_viirs_mat(Mat &img, Mat &lat, Mat &lon, bool sortoutput);
void resample_viirs(float **imgarr, float **latarr, float **lonarr, int nx, int ny, bool sortoutput);
void getsortingind(Mat &sind, int height);
void getadjustedsortingind(Mat &sind, const Mat &lat);
Mat resample_sort(const Mat &sind, const Mat &img);

// utils.cc
void	eprintf(const char *fmt, ...);
void dumpmat(const char *filename, Mat &m);
void dumpfloat(const char *filename, float *buf, int nbuf);
