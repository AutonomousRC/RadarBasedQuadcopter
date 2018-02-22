#ifndef __LFM_FFT_H__
#define __LFM_FFT_H__

#include <time.h>

#include <GL/glut.h>
#include <GL/gl.h>

void set_hanning(double *, int);
void draw_hanning_windowed_fft(void);
void hanning_window(void);

#endif
