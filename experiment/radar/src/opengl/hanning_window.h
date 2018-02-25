#ifndef __HANNING_WINDOW_H__
#define __HANNING_WINDOW_H__

#include <time.h>

#include <GL/glut.h>
#include <GL/gl.h>

void set_hanning(double *, int);
void draw_hanning_windowed_fft(void);
void hanning_window(void);

#endif
