#ifndef __LFM_FFT_H__
#define __LFM_FFT_H__

#include <time.h>

#include <GL/glut.h>
#include <GL/gl.h>

typedef GLfloat GLTmatrix[16];
typedef GLfloat GLTvector3[3];

#define GLT_PI          3.14159265358979323846
#define GLT_PI_DIV_180  0.017453292519943296
#define DEG2RAD(x)  ((x) * GLT_PI_DIV_180)

void glt_scale_vector(GLTvector3, const GLfloat);
GLfloat glt_get_vector_length_sqrt(const GLTvector3);
GLfloat glt_get_vector_length(const GLTvector3);
void glt_normalize_vector(GLTvector3);
void glt_draw_torus(GLfloat, GLfloat, GLint, GLint);
void glt_load_identity_matrix(GLTmatrix);
void glt_rotation_matrix(float, float, float, float, GLTmatrix);
void draw_torus(void);
void torus_reshape(int, int);
void torus_timer(int);

#endif
