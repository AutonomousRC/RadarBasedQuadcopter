#ifndef __HAMMING_WINDOW_H__
#define __HAMMING_WINDOW_H__

#include <time.h>

#include <GL/glut.h>
#include <GL/gl.h>

void set_hamming(double *, int);
void draw_hamming_windowed_fft(void);
void hamming_window(void);

#endif
