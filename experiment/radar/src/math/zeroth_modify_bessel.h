#ifndef __ZEROTH_MODIFY_BESSEL_H__
#define __ZEROTH_MODIFY_BESSEL_H__

#include "complex_math.h"
#include <stdbool.h>

#define NumBitsPerChar		8U

typedef struct {
	struct {
		unsigned int wordH;
		unsigned int wordL;
	} words;
} BigEndianIEEEDouble;

typedef struct {
	struct {
		unsigned int wordL;
		unsigned int wordH;
	} words;
} LittleEndianIEEEDouble;

typedef struct {
	union {
		float wordLreal;
		unsigned int wordLuint;
	} wordL;
} IEEESingle;

void zeroth_modify_bessel(double *, complex *, int);
void cbesi(complex, double, int, complex *, int *, int *);
int cseri(complex, double, int, int, complex *, double, double, double);
void gammaln(double *);
int cmlri(complex, double, int, int, complex *, double);
int casyi(complex, double, int, int, complex *, double, double, double);
int b_cuoik(complex, double, int, int, int, complex y[2], double, double, double);
int cuoik(complex, double, int, int, int, complex *, double, double, double);
void b_cuni2(complex, double, int, int, complex y[2], double, double, double, double, int *, int *);
void cuni1(complex, double, int, int, complex *, double, double, double, int *, int *);
void cuni2(complex, double, int, int, complex *, double, double, double, int *, int *);
void cbuni(complex, double, int, int, complex *, int, double, double, double, double, int *, int *);
int cwrsk(complex, double, int, int, complex *, complex cw[2], double, double, double);
int cuchk(complex, double, double);
void b_cunik(complex, double, int, int, double, int *, complex cwrk[16], complex *, complex *, complex *, complex *);
void cunik(complex, double, int, int, double, int, complex *, complex *, complex *);
void b_cunhj(complex, double, int, double, complex *, complex *, complex *, complex *, complex *, complex *);
void cunhj(complex, double, double, complex *, complex *, complex *, complex *);
int cacai(complex, double, int, int, complex *, double, double, double, double);
complex cairy(complex, int, int);
int cuchk(complex, double, double);
void b_log(complex *);
int b_ckscl(complex, int, complex y[2], double, double, double);
complex calccs(complex, complex, double);
int ckscl(complex, complex *, double);
int b_cbknu(complex, double, int, int, complex y[2], double, double, double);
void cbknu(complex, double, int, double, complex *, int *);
void b_cosh(complex *);
void b_sinh(complex *);
void b_exp(complex *);
void b_fix(double *);
double rt_hypotd_snf(double, double);
bool rtIsNaN(double);
double rt_atan2d_snf(double, double);
bool rtIsInf(double);
void b_sqrt(complex *);
float rtGetInfF(void);
float rtGetMinusInfF(void);
float rtGetNaNF(void);

#endif
