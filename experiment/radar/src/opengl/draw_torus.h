#ifndef __LFM_FFT_H__
#define __LFM_FFT_H__

#include <time.h>

#include <GL/glut.h>
#include <GL/gl.h>

typedef GLfloat GLTmatrix[16];
typedef GLfloat GLTvector3[3];

#define GLT_PI          3.14159265358979323846
#define GLT_PI_DIV_180  0.017453292519943296
#define gltDegToRad(x)  ((x) * GLT_PI_DIV_180)

void gltScaleVector(GLTvector3, const GLfloat);
GLfloat gltGetVectorLengthSqrd(const GLTvector3);
GLfloat gltGetVectorLength(const GLTvector3);
void gltNormalizeVector(GLTvector3);
void gltDrawTorus(GLfloat, GLfloat, GLint, GLint);
void gltLoadIdentityMatrix(GLTmatrix);
void gltRotationMatrix(float, float, float, float, GLTmatrix);
void draw_torus(void);

#endif
