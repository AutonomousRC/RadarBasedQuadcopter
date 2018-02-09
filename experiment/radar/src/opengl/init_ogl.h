#ifndef __INIT_OGL_H__
#define __INIT_OGL_H__

#include <time.h>

#include <GL/glut.h>
#include <GL/gl.h>

void ogl_init(int, char **, char *);
void ogl_post_process(void (*)(void));

#endif
