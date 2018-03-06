#ifndef __BLUESTEIN_FFT_H__
#define __BLUESTEIN_FFT_H__

#include "complex_math.h"

int calc_align_idx(int);
void bluestein_first(comp *, int, double *, double *, comp *, int);
void bluestein_second(comp *, double *, double *, comp *, int);
void bluestein_third(comp *, double *, double *, comp *, int);
void set_bluestein(comp *, int);
void bluestein_fft(comp *, comp *, int);
void fft_shift(double *);

#endif
