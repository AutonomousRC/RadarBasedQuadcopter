#ifndef __OGL_HELPER_H__
#define __OGL_HELPER_H__

#include <time.h>

#include <GL/glut.h>
#include <GL/gl.h>

void drawString(char *);
void drawStringBig(char *);
void slc(float, float, int, float *);
void init_screen(void);
void rect_screen(void);
#if 0
void p_arr(double *, int);
#endif
void display(void);
void reshape(int w, int h);
void draw_function(void);

#endif
