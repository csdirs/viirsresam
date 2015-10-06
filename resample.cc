//
// Resampling of image based on latitude
//

#include "viirsresam.h"
#include "sort.h"

// Pixels in deletion zone are given this value.
// It's -999.0 for band 13, and scaled value of integer 65533 for other bands.
// Set during entery point of resampling, in resample_viirs.
float DELETION_ZONE_VALUE = -999.0;

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

void
adjustbreakpoints(const Mat &slat, Mat &breakpointsT)
{
	CHECKMAT(slat, CV_32FC1);
	
	int nscans = slat.rows/NDETECTORS;
	short breakpoints[1+NCOLUMN_BREAKS] = {0, 5, 87, 170, 358, 567, 720, 850, 997, 1120, 1275, 1600};
	short detectorT[NCOLUMN_BREAKS-1] = {2, 8, 1, 2, 1, 2, 1, 2, 1, 0};
	
	breakpointsT = Mat::zeros(nscans, NCOLUMN_BREAKS, CV_32SC1);
	
	// copy original break points because we don't change all the break points
	for(int y = 0; y < nscans; y++){
		for(int x = 0; x < NCOLUMN_BREAKS; x++){
			breakpointsT.at<int>(y, x) = SORT_BREAK_POINTS[x];
		}
	}
	
	// Break points for 2nd scan to last scan.
	// N.B. Terminating break point (1600) is not set here.
	for(int j = 0; j < NCOLUMN_BREAKS-1; j++){
		int d = detectorT[j];
		int br = breakpoints[j+1];

		for(int k = 1; k < nscans; k++){
			const float *currow = slat.ptr<float>(k*NDETECTORS+d-1);
			const float *nextrow = slat.ptr<float>(k*NDETECTORS+d);

			int leftsign = SIGN(nextrow[br-1] - currow[br-1]);
			int rightsign = SIGN(nextrow[br+1] - currow[br+1]);

			// find if order is ascending (+1) or descending (-1)
			int order = SIGN(nextrow[slat.cols/2] - currow[slat.cols/2]);
			if(order == 0){
				continue;
			}
			
			if(leftsign == order && rightsign == order){
				breakpointsT.at<int>(k, j) = br;
				continue;
			}

			int signshift = -1;
			int signoff = leftsign;
			if(rightsign != order){
				signshift = +1;
				signoff = rightsign;
			}
			
			int count = 1;
			while(signoff != order && SIGN(breakpoints[j+1+signshift] - br+signshift*count) == signshift){
				count++;
				signoff = SIGN(nextrow[br+signshift*count] - currow[br+signshift*count]);
			}
			breakpointsT.at<int>(k, j) = br + signshift*(count-1);
		}
	}
}

// Generate a image of latitude sorting indices.
//
// sind -- sorting indices (output)
// height -- height of the output
//
void
getsortingind(Mat &sind, int height)
{
	sind = Mat::zeros(height, VIIRS_WIDTH, CV_32SC1);
	
	int x = 0;
	for(int i = 0; i < NCOLUMN_BREAKS; i++){
		int xe = SORT_BREAK_POINTS[i];
		for(; x < xe; x++){
			for(int y = 0; y < NDETECTORS; y++){
				sind.at<int>(y, x) = y + SORT_FIRST[y][i];
			}
			for(int y = NDETECTORS; y < height-NDETECTORS; y++){
				sind.at<int>(y, x) = y + SORT_MID[y%NDETECTORS][i];
			}
			for(int y = height-NDETECTORS; y < height; y++){
				sind.at<int>(y, x) = y + SORT_LAST[y%NDETECTORS][i];
			}
		}
	}
	
	x = VIIRS_WIDTH-1;
	for(int i = 0; i < NCOLUMN_BREAKS; i++){
		int xe = VIIRS_WIDTH - SORT_BREAK_POINTS[i];
		for(; x >= xe; x--){
			for(int y = 0; y < NDETECTORS; y++){
				sind.at<int>(y, x) = y + SORT_FIRST[y][i];
			}
			for(int y = NDETECTORS; y < height-NDETECTORS; y++){
				sind.at<int>(y, x) = y + SORT_MID[y%NDETECTORS][i];
			}
			for(int y = height-NDETECTORS; y < height; y++){
				sind.at<int>(y, x) = y + SORT_LAST[y%NDETECTORS][i];
			}
		}
	}
}

static inline void
getsortingind_left(Mat &sind, int startrow, int endrow, const Mat &breakpoints,
	short offset[][NCOLUMN_BREAKS])
{
	for(int y = startrow; y < endrow; y++){
		int x = 0;
		int scan = y/NDETECTORS;
		for(int i = 0; i < NCOLUMN_BREAKS; i++){
			int xe = breakpoints.at<int>(scan, i);
			for(; x < xe; x++){
				sind.at<int>(y, x) = y + offset[y%NDETECTORS][i];
			}
		}
	}
}

static inline void
getsortingind_right(Mat &sind, int startrow, int endrow, const Mat &breakpoints,
	short offset[][NCOLUMN_BREAKS])
{
	for(int y = startrow; y < endrow; y++){
		int x = VIIRS_WIDTH-1;
		int scan = y/NDETECTORS;
		for(int i = 0; i < NCOLUMN_BREAKS; i++){
			int xe = VIIRS_WIDTH - breakpoints.at<int>(scan, i);
			for(; x >= xe; x--){
				sind.at<int>(y, x) = y + offset[y%NDETECTORS][i];
			}
		}
	}
}

// Generate a image of latitude sorting indices with given breakpoints.
//
// sind -- sorting indices (output)
// height -- height of the output
// breakpoints -- breakpoints with shape [scan] x [number of breakpoints]+2
//
void
getsortingind1(Mat &sind, int height, const Mat &leftbreaks, const Mat &rightbreaks)
{
	CHECKMAT(leftbreaks, CV_32SC1);
	CHECKMAT(rightbreaks, CV_32SC1);
	CV_Assert(leftbreaks.cols == NCOLUMN_BREAKS);
	CV_Assert(rightbreaks.cols == NCOLUMN_BREAKS);

	getsortingind_left(sind, 0, NDETECTORS, leftbreaks, SORT_FIRST);
	getsortingind_left(sind, NDETECTORS, height-NDETECTORS, leftbreaks, SORT_MID);
	getsortingind_left(sind, height-NDETECTORS, height, leftbreaks, SORT_LAST);
	
	getsortingind_right(sind, 0, NDETECTORS, rightbreaks, SORT_FIRST);
	getsortingind_right(sind, NDETECTORS, height-NDETECTORS, rightbreaks, SORT_MID);
	getsortingind_right(sind, height-NDETECTORS, height, rightbreaks, SORT_LAST);
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
Mat
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
//
// T -- 3 SST values
// lat -- 3 latitudes
// lon -- 3 longitudes
// targlat -- target latitude
// targlon -- target longitude
// res -- spatial resolution
//
double
geoapprox(const float *T, const float *lat, const float *lon, float targlat, float targlon, double res)
{
	// none valid
	if(isinvalid(T[0]) && isinvalid(T[1]) && isinvalid(T[2]))
		return DELETION_ZONE_VALUE;

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

// Linear interpolation at x between points (x0, y0) and (x1, y1).
inline double
linearinterp(double x0, double y0, double x1, double y1, double x)
{
	double lam = (x - x0)/(x1 - x0);
	return (1-lam)*y0 + lam*y1;
}

// Interpolate longitude based on latitude sorting order.
// This makes the longitude monotonic.
//
// sind -- latitude sorting indices
// slon -- sorted longitude
// lon -- unsorted longitude
// n -- number of elements
// dst -- destination of interpolation (output)
//
void
interplon(const int *sind, const float *slon, const float *lon, int n, float *dst)
{
	vector<int> buf;
	int i;
	
	// extrapolate the reordered points before the first "kept order" point
	for(i = 0; i < n; i++){
		if(sind[i] == i){
			break;
		}
		buf.push_back(i);
	}
	for(int j = 0; j < (int)buf.size(); j++){
		dst[buf[j]] = slon[i];
	}
	buf.clear();
	double prevkeep = i;
	double prevlon = slon[i];
	
	// interpolate reordered points
	for(; i < n; i++){
		// sneak in middle of swath (between middle two detectors)
		if(i%NDETECTORS == NDETECTORS/2){
			double curkeep = i-0.5;
			double curlon = ((double)lon[i]+lon[i-1])/2.0;
			
			// interpolate at points in the buffer and clear the buffer
			for(int j = 0; j < (int)buf.size(); j++){
				int k = buf[j];
				dst[k] = linearinterp(prevkeep, prevlon, curkeep, curlon, k);
			}
			buf.clear();
			
			prevkeep = curkeep;
			prevlon = curlon;
		}
		
		if(sind[i] == i){	// kept order
			// interpolate at points in the buffer and clear the buffer
			for(int j = 0; j < (int)buf.size(); j++){
				int k = buf[j];
				dst[k] = linearinterp(prevkeep, prevlon, i, slon[i], k);
			}
			buf.clear();
			
			prevkeep = i;
			prevlon = slon[i];
			dst[i] = slon[i];
		}else{	// reordered
			buf.push_back(i);
		}
	}
	
	// extrapolate the reordered points after the last "kept order" point
	for(int j = 0; j < (int)buf.size(); j++){
		dst[buf[j]] = prevlon;
	}
}

// Resample 1D data.
//
// sind -- sorting indices
// sval -- sorted values
// slat -- sorted latitude
// slon -- sorted longitude
// ilon -- interpolated longitude
// n -- number of elements
// res -- spatial resolution
// rval -- resampled values (output)
//
static void
resample1d(const int *sind, const float *sval, const float *slat, const float *slon,
	const float *ilon, int n, double res, float *rval)
{
	int i;
	
	// Interpolate the middle values.
	// Set first and last values to an invalid value for now
	// because interpolation requires 3 consecutive values.
	rval[0] = INVALID_TEMP;
	for(i = 1; i < n-1; i++){
		//if(sind[i] == i){	// kept order
		//	rval[i] = sval[i];
		//}else{	// reordered
		//	rval[i] = geoapprox(&sval[i-1], &slat[i-1], &slon[i-1], slat[i], ilon[i], res);
		//}
		rval[i] = geoapprox(&sval[i-1], &slat[i-1], &slon[i-1], slat[i], ilon[i], res);
	}
	rval[n-1] = INVALID_TEMP;
	
	// extrapolate beginning values that are invalid
	for(i = 0; i < n; i++){
		if(!isinvalid(rval[i])){
			break;
		}
	}
	for(int k = 0; k < i; k++){
		rval[k] = rval[i];
	}
	
	// extrapolate ending values that are invalid
	for(i = n-1; i >= 0; i--){
		if(!isinvalid(rval[i])){
			break;
		}
	}
	for(int k = i+1; k < n; k++){
		rval[k] = rval[i];
	}
}

// Resample a 2D image.
//
// sortidx -- latitude sorting indices
// ssrc -- image to resample already sorted
// slat -- sorted latitude
// slon -- sorted longitude
// lon -- unsorted longitude
// dst -- resampled image (output)
// 
static void
resample2d(const Mat &sortidx, const Mat &ssrc, const Mat &slat, const Mat &slon,
	const Mat &lon, Mat &dst)
{
	Mat ilon;

	CHECKMAT(ssrc, CV_32FC1);
	CHECKMAT(slat, CV_32FC1);
	CHECKMAT(slon, CV_32FC1);
	CHECKMAT(sortidx, CV_32SC1);
	CV_Assert(ssrc.data != dst.data);

	int width = ssrc.cols;
	int height = ssrc.rows;
	
	// compute resolution per column based on the first two rows
	Mat _res = Mat::zeros(1, width, CV_64FC1);
	double *res = (double*)_res.data;
	//float *lat1 = (float*)slat.ptr(0);
	//float *lon1 = (float*)slon.ptr(0);
	//float *lat2 = (float*)slat.ptr(1);
	//float *lon2 = (float*)slon.ptr(1);
	for(int j = 0; j < width; j++){
		//res[j] = geodist(lat1[j], lon1[j], lat2[j], lon2[j])/4.0;
		double x = 2*j/(double)width - 1.0;
		res[j] = 0.2*SQ(x) + 0.2;
	}
	if(DEBUG)dumpmat("res.bin", _res);

	// allocate output and temporary bufferes for each column
	dst = Mat::zeros(height, width, CV_32FC1);
	if(DEBUG) ilon = Mat::zeros(height, width, CV_32FC1);
	Mat sindcol = Mat::zeros(height, 1, CV_32SC1);
	Mat ssrccol = Mat::zeros(height, 1, CV_32FC1);
	Mat slatcol = Mat::zeros(height, 1, CV_32FC1);
	Mat sloncol = Mat::zeros(height, 1, CV_32FC1);
	Mat loncol = Mat::zeros(height, 1, CV_32FC1);
	Mat dstcol = Mat::zeros(height, 1, CV_32FC1);
	Mat iloncol = Mat::zeros(height, 1, CV_32FC1);
	
	// resample each column
	for(int j = 0; j < width; j++){
		// copy columns to contiguous Mats, so we don't have to worry about stride
		sortidx.col(j).copyTo(sindcol.col(0));
		ssrc.col(j).copyTo(ssrccol.col(0));
		slat.col(j).copyTo(slatcol.col(0));
		slon.col(j).copyTo(sloncol.col(0));
		lon.col(j).copyTo(loncol.col(0));
		
		// interpolate longitude to make it monotonic
		interplon(sindcol.ptr<int>(0),
			sloncol.ptr<float>(0),
			loncol.ptr<float>(0),
			height,
			iloncol.ptr<float>(0));
		if(DEBUG) iloncol.col(0).copyTo(ilon.col(j));
		
		// resample and copy column to output
		resample1d(sindcol.ptr<int>(0),
			ssrccol.ptr<float>(0),
			slatcol.ptr<float>(0),
			sloncol.ptr<float>(0),
			iloncol.ptr<float>(0),
			height,
			res[j],
			dstcol.ptr<float>(0));
		dstcol.col(0).copyTo(dst.col(j));
	}
	if(DEBUG)dumpmat("ilon.bin", ilon);
}

void
resample_viirs_mat(Mat &img, Mat &lat, Mat &lon, float delval, bool sortoutput)
{
	Mat sind, dst;
	
	CHECKMAT(img, CV_32FC1);
	CHECKMAT(lat, CV_32FC1);
	CHECKMAT(lon, CV_32FC1);
	if(DEBUG)dumpmat("before.bin", img);
	if(DEBUG)dumpmat("lat.bin", lat);
	if(DEBUG)dumpmat("lon.bin", lon);
	
	if(DEBUG) printf("resampling debugging is turned on!\n");
	
	DELETION_ZONE_VALUE = delval;
	if(DEBUG) printf("deletion zone value is %f\n", DELETION_ZONE_VALUE);

	int ny = img.rows;
	int nx = img.cols;
	if(ny%NDETECTORS != 0){
		eprintf("invalid height %d (not multiple of %d)\n", ny, NDETECTORS);
	}
	if(nx != VIIRS_WIDTH){
		eprintf("invalid width %d; want %d", nx, VIIRS_WIDTH);
	}
	
	getsortingind(sind, ny);
	Mat slat = resample_sort(sind, lat);
	Mat slon = resample_sort(sind, lon);
	Mat simg = resample_sort(sind, img);
	if(DEBUG)dumpmat("sind.bin", sind);
	if(DEBUG)dumpmat("simg.bin", simg);
	if(DEBUG)dumpmat("slat.bin", slat);
	if(DEBUG)dumpmat("slon.bin", slon);
	
	resample2d(sind, simg, slat, slon, lon, dst);
	if(DEBUG)dumpmat("after.bin", dst);
	
	if(!sortoutput){
		dst = resample_unsort(sind, dst);
	}
	CV_Assert(dst.size() == img.size() && dst.type() == img.type());
	dst.copyTo(img);
	if(DEBUG)dumpmat("final.bin", img);
	
	if(sortoutput){
		slat.copyTo(lat);
		slon.copyTo(lon);
	}
	if(DEBUG)exit(3);
}

// Resample a VIIRS swath image.
//
// _img -- swath brightness temperature image to be resampled (input & output)
// _lat -- corresponding latitude image (input & output)
// _lon -- corresponding longitude image (input & output)
// nx -- width of image (should be 3200 for VIIRS)
// ny -- height of image (5408 or 5392 for ~10 min VIIRS granule)
// delval -- value used for deletion zone in _img
// sortoutput -- indicates if output should be in latitude sorted order
//
void
resample_viirs(float **_img, float **_lat, float **_lon, int nx, int ny, float delval, bool sortoutput)
{
	Mat _sind, sind, leftbreaks, rightbreaks, dst;

	if(DEBUG) printf("resampling debugging is turned on!\n");
	
	DELETION_ZONE_VALUE = delval;
	if(DEBUG) printf("deletion zone value is %f\n", DELETION_ZONE_VALUE);

	if(ny%NDETECTORS != 0){
		eprintf("invalid height %d (not multiple of %d)\n", ny, NDETECTORS);
	}
	if(nx != VIIRS_WIDTH){
		eprintf("invalid width %d; want %d", nx, VIIRS_WIDTH);
	}
	
	// Mat wrapper around external buffer.
	// Caller of this function still reponsible for freeing the buffers.
	Mat img(ny, nx, CV_32FC1, &_img[0][0]);
	Mat lat(ny, nx, CV_32FC1, &_lat[0][0]);
	Mat lon(ny, nx, CV_32FC1, &_lon[0][0]);
	if(DEBUG)dumpmat("before.bin", img);
	if(DEBUG)dumpmat("lat.bin", lat);
	if(DEBUG)dumpmat("lon.bin", lon);

	getsortingind(_sind, ny);
	if(DEBUG)dumpmat("_sind.bin", _sind);
	Mat _slat = resample_sort(_sind, lat);
	if(DEBUG)dumpmat("_slat.bin", _slat);
	if(true){
		sind = Mat::zeros(ny, VIIRS_WIDTH, CV_32SC1);
		
		// left half
		adjustbreakpoints(_slat, leftbreaks);
		
		// right half
		Mat slatflipped = Mat::zeros(_slat.size(), _slat.type());
		flip(_slat, slatflipped, 1);
		if(DEBUG)dumpmat("slatflipped.bin", slatflipped);
		adjustbreakpoints(slatflipped, rightbreaks);
		
		if(DEBUG)dumpmat("leftbreaks.bin", leftbreaks);
		if(DEBUG)dumpmat("rightbreaks.bin", rightbreaks);
		getsortingind1(sind, ny, leftbreaks, rightbreaks);
	}else{
		Mat testBP = Mat::zeros(_slat.rows/NDETECTORS, NCOLUMN_BREAKS, CV_32SC1);
		for(int i = 0; i < testBP.rows; i++){
			for(int j = 0; j < testBP.cols; j++){
				testBP.at<int>(i, j) = SORT_BREAK_POINTS[j];
			}
		}
		getsortingind1(sind, ny, testBP, testBP);
	}

	Mat slat = resample_sort(sind, lat);
	Mat slon = resample_sort(sind, lon);
	Mat simg = resample_sort(sind, img);
	if(DEBUG)dumpmat("sind.bin", sind);
	if(DEBUG)dumpmat("simg.bin", simg);
	if(DEBUG)dumpmat("slat.bin", slat);
	if(DEBUG)dumpmat("slon.bin", slon);
	
	resample2d(sind, simg, slat, slon, lon, dst);
	if(DEBUG)dumpmat("after.bin", dst);
	
	if(!sortoutput){
		dst = resample_unsort(sind, dst);
	}
	CV_Assert(dst.size() == img.size() && dst.type() == img.type());
	dst.copyTo(img);
	if(DEBUG)dumpfloat("final.bin", &_img[0][0], nx*ny);
	
	if(sortoutput){
		slat.copyTo(lat);
		slon.copyTo(lon);
	}
	if(DEBUG)exit(3);
}
