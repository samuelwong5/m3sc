#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define _STATIC_SF

#ifdef _STATIC_SF
void sf_init(int);
void sf_release(void);
#endif
double sf_index(int, int, double *);
double * SFactors(int);

int FastSN(double *, double *, double *, double *, int, int);
int FastTN(double *, double *, double *, double *, int, int);
int FastUN(double *, double *, double *, double *, int, int);

int SN(double *, double *, int, double *);
int TN(double *, double *, int, double *);
int UN(double *, double *, int, double *);