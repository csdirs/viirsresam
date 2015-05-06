// allocate_2d.cc
float ** allocate_2d_f(int n1, int n2);
int ** allocate_2d_i(int n1, int n2);

// readwrite_viirs.cc
int readwrite_viirs(unsigned short **buffer, unsigned long long * dimsizes, float * gain, float * offset, 
                    char * filename, char * BTstr, int readwrite);
int readwrite_viirs_float(float **buffer, unsigned long long * dimsizes, char * filename, char * BTstr, int readwrite);
int write_viirs_destriping_attribute(char * filename, char * attrFieldStr, char * attrNameStr, float destrval);

// resample_viirs.c
void resample_viirs(float **imgarr, float **latarr, float **lonarr, int nx, int ny, float min, float max);
