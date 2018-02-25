#ifndef __KAISER_WINDOW_H__
#define __KAISER_WINDOW_H__

#include <time.h>

#include <GL/glut.h>
#include <GL/gl.h>

void set_kaiser(double *, int);
void draw_kaiser_windowed_fft(void);
void kaiser_window(void);

#endif
