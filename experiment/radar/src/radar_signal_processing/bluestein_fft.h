#ifndef __BLUESTEIN_FFT_H__
#define __BLUESTEIN_FFT_H__

#include "complex_math.h"

void bluestein_first(complex *, int, double *, double *, complex *);
void bluestein_second(complex *, double *, double *, complex *);
void set_bluestein(complex *);
void fft(complex *, complex *);

#endif
