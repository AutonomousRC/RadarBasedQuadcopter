#ifndef __BLUESTEIN_FFT_H__
#define __BLUESTEIN_FFT_H__

#include "complex_math.h"

int calc_align_idx(int);
void bluestein_first(complex *, int, double *, double *, complex *, int);
void bluestein_second(complex *, double *, double *, complex *, int);
void bluestein_third(complex *, double *, double *, complex *, int);
void set_bluestein(complex *, int);
void bluestein_fft(complex *, complex *, int);
void fft_shift(double *);

#endif
