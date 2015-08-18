//
// Resampling of image based on latitude
//

#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <math.h>
#include "viirsresam.h"

#define SGN(A)   ((A) > 0 ? 1 : ((A) < 0 ? -1 : 0 ))
#define SQ(x)	((x)*(x))
#define CHECKMAT(M, T)	CV_Assert((M).type() == (T) && (M).isContinuous())

enum {
	VIIRS_SWATH_SIZE = 16,
	MAX_TEMP = 350,	// in Kelvin
	MIN_TEMP = 0,	// in Kelvin
	DEBUG = true,
};

using namespace cv;

inline bool
isinvalid(float x)
{
	return x > MAX_TEMP || x < MIN_TEMP;
}

static void
eprintf(const char *fmt, ...)
{
	va_list args;

	fflush(stdout);
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);

	if(fmt[0] != '\0' && fmt[strlen(fmt)-1] == ':')
		fprintf(stderr, " %s", strerror(errno));
	fprintf(stderr, "\n");

	exit(2);
}

static const char*
type2str(int type)
{
	switch(type) {
	default:
		return "UnknownType";
		break;
	case CV_8UC1:
		return "CV_8UC1";
		break;
	case CV_8SC1:
		return "CV_8SC1";
		break;
	case CV_16UC1:
		return "CV_16UC1";
		break;
	case CV_16SC1:
		return "CV_16SC1";
		break;
	case CV_32SC1:
		return "CV_32SC1";
		break;
	case CV_32FC1:
		return "CV_32FC1";
		break;
	case CV_64FC1:
		return "CV_64FC1";
		break;
	}
}

static void
dumpmat(const char *filename, Mat &m)
{
	int n;
	FILE *f;

	if(!m.isContinuous()) {
		eprintf("m not continuous");
	}
	f = fopen(filename, "w");
	if(!f) {
		eprintf("open %s failed:", filename);
	}
	n = fwrite(m.data, m.elemSize1(), m.rows*m.cols, f);
	if(n != m.rows*m.cols) {
		fclose(f);
		eprintf("wrote %d/%d items; write failed:", n, m.rows*m.cols);
	}
	fclose(f);
}

static void
dumpfloat(const char *filename, float *buf, int nbuf)
{
	int n;
	FILE *f;

	f = fopen(filename, "w");
	if(!f) {
		eprintf("open %s failed:", filename);
	}
	n = fwrite(buf, sizeof(*buf), nbuf, f);
	if(n != nbuf) {
		fclose(f);
		eprintf("wrote %d/%d items; write failed:", n, nbuf);
	}
	fclose(f);
}

template <class T>
static Mat
resample_unsort_(const Mat &sind, const Mat &img)
{
	Mat newimg;
	int i, j, k;
	int32_t *sp;
	T *ip;

	CHECKMAT(sind, CV_32SC1);
	CV_Assert(img.channels() == 1);

	newimg = Mat::zeros(img.rows, img.cols, img.type());
	sp = (int32_t*)sind.data;
	ip = (T*)img.data;
	k = 0;
	for(i = 0; i < newimg.rows; i++) {
		for(j = 0; j < newimg.cols; j++) {
			newimg.at<T>(sp[k], j) = ip[k];
			k++;
		}
	}
	return newimg;
}

// Returns the unsorted image of the sorted image img.
// Sind is the image of sort indices.
static Mat
resample_unsort(const Mat &sind, const Mat &img)
{
	switch(img.type()) {
	default:
		eprintf("unsupported type %s\n", type2str(img.type()));
		break;
	case CV_8UC1:
		return resample_unsort_<uchar>(sind, img);
		break;
	case CV_32FC1:
		return resample_unsort_<float>(sind, img);
		break;
	case CV_64FC1:
		return resample_unsort_<double>(sind, img);
		break;
	}
	// not reached
	return Mat();
}

template <class T>
static Mat
resample_sort_(const Mat &sind, const Mat &img)
{
	Mat newimg;
	int i, j, k;
	int32_t *sp;
	T *np;

	CHECKMAT(sind, CV_32SC1);
	CV_Assert(img.channels() == 1);

	newimg = Mat::zeros(img.rows, img.cols, img.type());
	sp = (int*)sind.data;
	np = (T*)newimg.data;
	k = 0;
	for(i = 0; i < newimg.rows; i++) {
		for(j = 0; j < newimg.cols; j++) {
			np[k] = img.at<T>(sp[k], j);
			k++;
		}
	}
	return newimg;
}

// Returns the sorted image of the unsorted image img.
// Sind is the image of sort indices.
static Mat
resample_sort(const Mat &sind, const Mat &img)
{
	switch(img.type()) {
	default:
		eprintf("unsupported type %s\n", type2str(img.type()));
		break;
	case CV_8UC1:
		return resample_sort_<uchar>(sind, img);
		break;
	case CV_32FC1:
		return resample_sort_<float>(sind, img);
		break;
	case CV_64FC1:
		return resample_sort_<double>(sind, img);
		break;
	}
	// not reached
	return Mat();
}


enum Pole {
	NORTHPOLE,
	SOUTHPOLE,
	NOPOLE,
};
typedef enum Pole Pole;

// Argsort latitude image 'lat' with given swath size.
// Image of sort indices are return in 'sortidx'.
static void
argsortlat(const Mat &lat, int swathsize, Mat &sortidx)
{
	int i, j, off, width, height, dir, d, split;
	Pole pole;
	Mat col, idx, botidx;
	Range colrg, toprg, botrg;

	CHECKMAT(lat, CV_32FC1);
	CV_Assert(swathsize >= 2);
	CV_Assert(lat.data != sortidx.data);

	width = lat.cols;
	height = lat.rows;
	sortidx.create(height, width, CV_32SC1);

	// For a column in latitude image, look at every 'swathsize' pixels
	// starting from 'off'. If they increases and then decreases, or
	// decreases and then increases, we're at the polar region.
	off = swathsize/2;

	pole = NOPOLE;

	for(j = 0; j < width; j++) {
		col = lat.col(j);

		// find initial direction -- increase, decrease or no change
		dir = 0;
		for(i = off+swathsize; i < height; i += swathsize) {
			dir = SGN(col.at<float>(i) - col.at<float>(i-swathsize));
			if(dir != 0)
				break;
		}

		// find change in direction if there is one
		for(; i < height; i += swathsize) {
			d = SGN(col.at<float>(i) - col.at<float>(i-swathsize));
			if(dir == 1 && d == -1) {
				CV_Assert(pole == NOPOLE || pole == NORTHPOLE);
				pole = NORTHPOLE;
				break;
			}
			if(dir == -1 && d == 1) {
				CV_Assert(pole == NOPOLE || pole == SOUTHPOLE);
				pole = SOUTHPOLE;
				break;
			}
		}

		if(i >= height) {
			pole = NOPOLE;
			if(dir >= 0)
				sortIdx(col, sortidx.col(j), CV_SORT_EVERY_COLUMN + CV_SORT_ASCENDING);
			else
				sortIdx(col, sortidx.col(j), CV_SORT_EVERY_COLUMN + CV_SORT_DESCENDING);
			continue;
		}

		split = i-swathsize;	// split before change in direction
		colrg = Range(j, j+1);
		toprg = Range(0, split);
		botrg = Range(split, height);

		if(pole == NORTHPOLE) {
			botidx = sortidx(botrg, colrg);
			sortIdx(col.rowRange(toprg), sortidx(toprg, colrg),
			        CV_SORT_EVERY_COLUMN + CV_SORT_ASCENDING);
			sortIdx(col.rowRange(botrg), botidx,
			        CV_SORT_EVERY_COLUMN + CV_SORT_DESCENDING);
			botidx += split;
		} else {	// pole == SOUTHPOLE
			botidx = sortidx(botrg, colrg);
			sortIdx(col.rowRange(toprg), sortidx(toprg, colrg),
			        CV_SORT_EVERY_COLUMN + CV_SORT_DESCENDING);
			sortIdx(col.rowRange(botrg), botidx,
			        CV_SORT_EVERY_COLUMN + CV_SORT_ASCENDING);
			botidx += split;
		}
	}
}

// Find distance between (lat1, lon1) and (lat2, lon2).
// Distances are computed using haversine formula 
// Other versions:
//	http://en.wikipedia.org/wiki/Great-circle_distance
//	http://www.movable-type.co.uk/scripts/latlong.html
//
double
geodist(double lat1, double lon1, double lat2, double lon2)
{
	const double R = 6371.0;
	double phi1 = (M_PI * lat1) / 180.0;
	double phi2 = (M_PI * lat2) / 180.0;
	double lam1 = (M_PI * lon1) / 180.0;
	double lam2 = (M_PI * lon2) / 180.0;
	double delta_phi = phi1 - phi2;
	double delta_lam = lam1 - lam2;
	
	return R*sqrt(SQ(cos((phi1+phi2)/2) * delta_lam) + SQ(delta_phi));
}

// Approximate from three (possibly invalid) values at lat/lon pairs.
double
geoapprox(const float *T, const float *lat, const float *lon, float targlat, float targlon, double res)
{
	// none valid
	if(isinvalid(T[0]) && isinvalid(T[1]) && isinvalid(T[2]))
		return -999;

	// one valid
	if(isinvalid(T[0]) && isinvalid(T[1]))
		return T[2];
	if(isinvalid(T[0]) && isinvalid(T[2]))
		return T[1];
	if(isinvalid(T[1]) && isinvalid(T[2]))
		return T[0];
	
	// at least two valid
	double sqres = SQ(res);
	double num = 0;
	double denom = 0;
	for(int i = 0; i < 3; i++){
		if(!isinvalid(T[i])){
			double d = geodist(targlat, targlon, lat[i], lon[i]);
			double w = exp(-SQ(d) / sqres);
			num += T[i] * w;
			denom += w;
		}
	}
	return num/denom;
}

// Resample 1D data.
//
// sind -- sorting indices
// slat -- sorted latitude
// sval -- sorted values
// dir -- diff direction (-1 or 1)
// rval -- resampled values (intput & output)
//
static void
resample1d(const int *sind, const float *slat, const float *slon, const float *sval, int n, double res, float *rval)
{
	int i;
	
	// copy first non-nan value for first row
	for(i = 0; i < n-1; i += 1){
		if(!isinvalid(sval[i])){
			rval[0] = sval[i];
			break;
		}
	}
	
	// interpolate the middle values
	for(i = 1; i < n-1; i += 1){
		if(sind[i] == i){	// kept order
			rval[i] = sval[i];
		}else{	// reordered
			// TODO: interpolate lon
			rval[i] = geoapprox(&sval[i-1], &slat[i-1], &slon[i-1], slat[i], slon[i], res);
		}
	}
	
	// copy last non-nan value to last row
	for(int k = i; k >= 0; k -= 1){
		if(!isinvalid(sval[k])){
			rval[i] = sval[k];
			break;
		}
	}
}

// Resample a 2D image.
//
// ssrc -- image to resample already sorted
// slat -- sorted latitude
// sortidx -- lat sorting indices
// dst -- resampled image (output)
// 
static void
resample2d(const Mat &ssrc, const Mat &slat, const Mat &slon, const Mat &sortidx, Mat &dst)
{
	int width, height;
	Mat col, idx, botidx;
	Range colrg, toprg, botrg;

	CHECKMAT(ssrc, CV_32FC1);
	CHECKMAT(slat, CV_32FC1);
	CHECKMAT(slon, CV_32FC1);
	CHECKMAT(sortidx, CV_32SC1);
	CV_Assert(ssrc.data != dst.data);

	width = ssrc.cols;
	height = ssrc.rows;
	
	// compute resolution per column based on the first two rows
	Mat _res = Mat::zeros(1, width, CV_64FC1);
	double *res = (double*)_res.data;
	float *lat1 = (float*)slat.ptr(0);
	float *lon1 = (float*)slon.ptr(0);
	float *lat2 = (float*)slat.ptr(1);
	float *lon2 = (float*)slon.ptr(1);
	for(int j = 0; j < width; j++){
		res[j] = geodist(lat1[j], lon1[j], lat2[j], lon2[j]);
	}
	if(DEBUG)dumpmat("res.bin", _res);

	dst = Mat::zeros(height, width, CV_32FC1);	// resampled values
	Mat sindcol = Mat::zeros(height, 1, CV_32SC1);
	Mat ssrccol = Mat::zeros(height, 1, CV_32FC1);
	Mat slatcol = Mat::zeros(height, 1, CV_32FC1);
	Mat sloncol = Mat::zeros(height, 1, CV_32FC1);
	Mat dstcol = Mat::zeros(height, 1, CV_32FC1);
	
	// resample each column
	for(int j = 0; j < width; j++){
		// copy columns to contiguous Mats, so we don't have to worry about stride
		sortidx.col(j).copyTo(sindcol.col(0));
		slat.col(j).copyTo(slatcol.col(0));
		slon.col(j).copyTo(sloncol.col(0));
		ssrc.col(j).copyTo(ssrccol.col(0));
		
		resample1d(sindcol.ptr<int>(0),
			slatcol.ptr<float>(0),
			sloncol.ptr<float>(0),
			ssrccol.ptr<float>(0),
			height,
			res[j],
			dstcol.ptr<float>(0));
		
		// copy resampled column to destination
		dstcol.col(0).copyTo(dst.col(j));
	}
}

// Resample VIIRS swatch image _img with corresponding
// latitude image _lat.
// _img[0..ny][0..nx]  - original image (brightness temperature)
// _lat[0..ny][0..nx]  - original latitude
// _lon[0..ny][0..nx]  - original longitude
// nx = width of image (should be 3200 for VIIRS)
// ny = height of image ( 5408 or 5392 for ~10 min VIIRS granule)
//
// TODO: use min, max arguments
void
resample_viirs(float **_img, float **_lat, float **_lon, int nx, int ny, float min, float max)
{
	Mat sind, dst;

	if(DEBUG) printf("resampling debugging is turned on!\n");
	
	// Mat wrapper around external buffer.
	// Caller of this function still reponsible for freeing the buffers.
	Mat img(ny, nx, CV_32FC1, &_img[0][0]);
	Mat lat(ny, nx, CV_32FC1, &_lat[0][0]);
	Mat lon(ny, nx, CV_32FC1, &_lon[0][0]);
	if(DEBUG)dumpmat("before.bin", img);
	if(DEBUG)dumpmat("lat.bin", lat);

	argsortlat(lat, VIIRS_SWATH_SIZE, sind);
	Mat slat = resample_sort(sind, lat);
	Mat slon = resample_sort(sind, lon);
	Mat simg = resample_sort(sind, img);
	if(DEBUG)dumpmat("sind.bin", sind);
	if(DEBUG)dumpmat("simg.bin", simg);
	if(DEBUG)dumpmat("slat.bin", slat);
	
	resample2d(simg, slat, slon, sind, dst);
	if(DEBUG)dumpmat("after.bin", dst);
	
	//simg = resample_interp(simg, lat);
	//if(DEBUG)dumpmat("simg3.bin", simg);
	//simg = resample_unsort(sind, simg);
	//if(DEBUG)dumpmat("simg4.bin", simg);

	//CV_Assert(simg.size() == img.size() && simg.type() == img.type());
	//simg.copyTo(img);
	//if(DEBUG)dumpfloat("final.bin", &_img[0][0], nx*ny);
	if(DEBUG)exit(3);
}
