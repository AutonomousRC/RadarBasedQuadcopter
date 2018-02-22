#ifndef __LFM_FFT_H__
#define __LFM_FFT_H__

#include <time.h>

#include <GL/glut.h>
#include <GL/gl.h>

void set_hamming(double *, int);
void draw_hamming_windowed_fft(void);
void hamming_window(void);

#endif
