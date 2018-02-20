#ifndef __INIT_OGL_H__
#define __INIT_OGL_H__

#include <time.h>

#include <GL/glut.h>
#include <GL/gl.h>

void ogl_init(int, char **, char *);
void ogl_3d_init(int, char **, char *);
void setup_rc(void);
void ogl_post_process(void (*)(void));
void ogl_3d_post_process(void (*)(void), void (*)(int, int), void (*)(int));

#endif
