#include <stdlib.h>


// these are all routines for allocation of a 2d array arr of size arr[n1][n2] for different data types
// of course, C++ template would be more elegant, but this way we can get by with just C compiler

char ** allocate_2d_c(int n1, int n2)
{
	int i;
	char ** p = NULL;
	p = (char **) malloc(n1*sizeof(char *));
	if(p==NULL) return NULL;
	p[0] = (char *) malloc(n1*n2*sizeof(char));
	if(p[0]==NULL) return NULL;
	for(i=1; i<n1; i++) {
		p[i] = &(p[0][i*n2]);
	}
	return p;
};

float ** allocate_2d_f(int n1, int n2)
{
	int i;
	float ** p = NULL;
	p = (float **) malloc(n1*sizeof(float *));
	if(p==NULL) return NULL;
	p[0] = (float *) malloc(n1*n2*sizeof(float));
	if(p[0]==NULL) return NULL;
	for(i=1; i<n1; i++) {
		p[i] = &(p[0][i*n2]);
	}
	return p;
};

double ** allocate_2d_d(int n1, int n2)
{
	int i;
	double ** p = NULL;
	p = (double **) malloc(n1*sizeof(double *));
	if(p==NULL) return NULL;
	p[0] = (double *) malloc(n1*n2*sizeof(double));
	if(p[0]==NULL) return NULL;
	for(i=1; i<n1; i++) {
		p[i] = &(p[0][i*n2]);
	}
	return p;
};

int ** allocate_2d_i(int n1, int n2)
{
	int i;
	int ** p = NULL;
	p = (int **) malloc(n1*sizeof(int *));
	if(p==NULL) return NULL;
	p[0] = (int *) malloc(n1*n2*sizeof(int));
	if(p[0]==NULL) return NULL;
	for(i=1; i<n1; i++) {
		p[i] = &(p[0][i*n2]);
	}
	return p;
};

short ** allocate_2d_s(int n1, int n2)
{
	int i;
	short ** p = NULL;
	p = (short **) malloc(n1*sizeof(short *));
	if(p==NULL) return NULL;
	p[0] = (short *) malloc(n1*n2*sizeof(short));
	if(p[0]==NULL) return NULL;
	for(i=1; i<n1; i++) {
		p[i] = &(p[0][i*n2]);
	}
	return p;
};

unsigned short ** allocate_2d_us(int n1, int n2)
{
	int i;
	unsigned short ** p = NULL;
	p = (unsigned short **) malloc(n1*sizeof(unsigned short *));
	if(p==NULL) return NULL;
	p[0] = (unsigned short *) malloc(n1*n2*sizeof(unsigned short));
	if(p[0]==NULL) return NULL;
	for(i=1; i<n1; i++) {
		p[i] = &(p[0][i*n2]);
	}
	return p;
};
