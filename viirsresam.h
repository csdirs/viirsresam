#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <math.h>
#include <opencv2/opencv.hpp>

using namespace cv;

#define nelem(x)	(sizeof(x)/sizeof((x)[0]))
#define CHECKMAT(M, T)	CV_Assert((M).type() == (T) && (M).isContinuous())
#define SQ(x)	((x)*(x))

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
int readwrite_viirs_float(float **buffer, unsigned long long * dimsizes, char * filename, char * BTstr, int readwrite);
int write_viirs_destriping_attribute(char * filename, char * attrFieldStr, char * attrNameStr, float destrval);

// resample.cc
void resample_viirs(float **imgarr, float **latarr, float **lonarr, int nx, int ny, float delval, bool sortoutput);
