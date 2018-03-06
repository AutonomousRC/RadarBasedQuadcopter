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

void zeroth_modify_bessel(double *, comp *, int);
void cbesi(comp, double, int, comp *, int *, int *);
int cseri(comp, double, int, int, comp *, double, double, double);
void gammaln(double *);
int cmlri(comp, double, int, int, comp *, double);
int casyi(comp, double, int, int, comp *, double, double, double);
int b_cuoik(comp, double, int, int, int, comp y[2], double, double, double);
int cuoik(comp, double, int, int, int, comp *, double, double, double);
void b_cuni2(comp, double, int, int, comp y[2], double, double, double, double, int *, int *);
void cuni1(comp, double, int, int, comp *, double, double, double, int *, int *);
void cuni2(comp, double, int, int, comp *, double, double, double, int *, int *);
void cbuni(comp, double, int, int, comp *, int, double, double, double, double, int *, int *);
int cwrsk(comp, double, int, int, comp *, comp cw[2], double, double, double);
int cuchk(comp, double, double);
void b_cunik(comp, double, int, int, double, int *, comp cwrk[16], comp *, comp *, comp *, comp *);
void cunik(comp, double, int, int, double, int, comp *, comp *, comp *);
void b_cunhj(comp, double, int, double, comp *, comp *, comp *, comp *, comp *, comp *);
void cunhj(comp, double, double, comp *, comp *, comp *, comp *);
int cacai(comp, double, int, int, comp *, double, double, double, double);
comp cairy(comp, int, int);
int cuchk(comp, double, double);
void b_log(comp *);
int b_ckscl(comp, int, comp y[2], double, double, double);
comp calccs(comp, comp, double);
int ckscl(comp, comp *, double);
int b_cbknu(comp, double, int, int, comp y[2], double, double, double);
void cbknu(comp, double, int, double, comp *, int *);
void b_cosh(comp *);
void b_sinh(comp *);
void b_exp(comp *);
void b_fix(double *);
double rt_hypotd_snf(double, double);
bool rtIsNaN(double);
double rt_atan2d_snf(double, double);
bool rtIsInf(double);
void b_sqrt(comp *);
float rtGetInfF(void);
float rtGetMinusInfF(void);
float rtGetNaNF(void);

#endif
