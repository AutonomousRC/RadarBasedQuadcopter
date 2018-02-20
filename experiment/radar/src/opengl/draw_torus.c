#include "init_ogl.h"
#include "ogl_helper.h"
#include "draw_torus.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if 0
typedef GLfloat GLTmatrix[16];
typedef GLfloat GLTvector3[3];

#define GLT_PI		3.14159265358979323846
#define GLT_PI_DIV_180	0.017453292519943296
#define gltDegToRad(x)	((x) * GLT_PI_DIV_180)
#endif

void gltScaleVector(GLTvector3 vVector, const GLfloat fScale)
{
	vVector[0] *= fScale;
	vVector[1] *= fScale;
	vVector[2] *= fScale;
}

GLfloat gltGetVectorLengthSqrd(const GLTvector3 vVector)
{
	return pow(vVector[0], 2) + pow(vVector[1], 2) + pow(vVector[2], 2);
}

GLfloat gltGetVectorLength(const GLTvector3 vVector)
{
	return (GLfloat)sqrt(gltGetVectorLengthSqrd(vVector));
}

void gltNormalizeVector(GLTvector3 vNormal)
{
	GLfloat fLength = 1.0f / gltGetVectorLength(vNormal);
	gltScaleVector(vNormal, fLength);
}

void gltDrawTorus(GLfloat majorRadius, GLfloat minorRadius, GLint numMajor, GLint numMinor)
{
	GLTvector3 vNormal;
	double majorStep = 2.0f * GLT_PI / numMajor;
	double minorStep = 2.0f * GLT_PI / numMinor;
	int i, j;

	for(i = 0; i < numMajor; ++i)
	{
		double a0 = i * majorStep;
		double a1 = a0 + majorStep;

		GLfloat x0 = (GLfloat)cos(a0);
		GLfloat y0 = (GLfloat)sin(a0);
		GLfloat x1 = (GLfloat)cos(a1);
		GLfloat y1 = (GLfloat)sin(a1);

		glBegin(GL_TRIANGLE_STRIP);

		for(j = 0; j <= numMinor; ++j)
		{
			double b = j * minorStep;
			GLfloat c = (GLfloat)cos(b);
			GLfloat r = minorRadius * c + majorRadius;
			GLfloat z = minorRadius * (GLfloat)sin(b);
			glTexCoord2f((float)(i) / (float)(numMajor), (float)(j) / (float)(numMinor));

			vNormal[0] = x0 * c;
			vNormal[1] = y0 * c;
			vNormal[2] = z / minorRadius;

			glVertex3f(x0 * r, y0 * r, z);
			glTexCoord2f((float)(i + 1) / (float)(numMajor), (float)(j) / (float)(numMinor));

			vNormal[0] = x1 * c;
			vNormal[1] = y1 * c;
			vNormal[2] = z / minorRadius;

			glNormal3fv(vNormal);
			glVertex3f(x1 * r, y1 * r, z);
		}

		glEnd();
	}
}

void gltLoadIdentityMatrix(GLTmatrix m)
{
	static GLTmatrix identity = {	1.0f, 0.0f, 0.0f, 0.0f,
					0.0f, 1.0f, 0.0f, 0.0f,
					0.0f, 0.0f, 1.0f, 0.0f,
					0.0f, 0.0f, 0.0f, 1.0f };
	memcpy(m, identity, sizeof(GLTmatrix));
}

void gltRotationMatrix(float angle, float x, float y, float z, GLTmatrix mMatrix)
{
	float vecLength, sinSave, cosSave, oneMinusCos;
	float xx, yy, zz, xy, yz, zx, xs, ys, zs;

	if(x == 0.0f && y == 0.0f && z == 0.0f)
	{
		gltLoadIdentityMatrix(mMatrix);
		return;
	}

	vecLength = (float)sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

	x /= vecLength;
	y /= vecLength;
	z /= vecLength;

	sinSave = (float)sin(angle);
	cosSave = (float)cos(angle);
	oneMinusCos = 1.0f - cosSave;

	xx = pow(x, 2);
	yy = pow(y, 2);
	zz = pow(z, 2);
	xy = x * y;
	yz = y * z;
	zx = z * x;
	xs = x * sinSave;
	ys = y * sinSave;
	zs = z * sinSave;

	mMatrix[0] = (oneMinusCos * xx) + cosSave;
	mMatrix[4] = (oneMinusCos * xy) - zs;
	mMatrix[8] = (oneMinusCos * zx) + ys;
	mMatrix[12] = 0.0f;

	mMatrix[1] = (oneMinusCos * xy) + zs;
	mMatrix[5] = (oneMinusCos * yy) + cosSave;
	mMatrix[9] = (oneMinusCos * yz) - xs;
	mMatrix[13] = 0.0f;

	mMatrix[2] = (oneMinusCos * zx) - ys;
	mMatrix[6] = (oneMinusCos * yz) + xs;
	mMatrix[10] = (oneMinusCos * zz) + cosSave;
	mMatrix[14] = 0.0f;

	mMatrix[3] = 0.0f;
	mMatrix[7] = 0.0f;
	mMatrix[11] = 0.0f;
	mMatrix[15] = 1.0f;
}

void draw_torus(void)
{
	GLTmatrix rot_mat, translation_mat, transform_mat;

	static GLfloat y_rot = 0.0f;
	y_rot += 0.5f;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	gltRotationMatrix(gltDegToRad(y_rot), 0.0f, 1.0f, 0.f, transform_mat);
	transform_mat[12] = 0.0f;
	transform_mat[13] = 0.0f;
	transform_mat[14] = -2.5f;

	glLoadMatrixf(transform_mat);

	gltDrawTorus(0.35f, 0.15f, 40, 20);

	glutSwapBuffers();
}
