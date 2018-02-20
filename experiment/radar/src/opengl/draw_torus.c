#include "init_ogl.h"
#include "ogl_helper.h"
#include "draw_torus.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void glt_scale_vector(GLTvector3 v_vector, const GLfloat f_scale)
{
	v_vector[0] *= f_scale;
	v_vector[1] *= f_scale;
	v_vector[2] *= f_scale;
}

GLfloat glt_get_vector_length_sqrt(const GLTvector3 v_vector)
{
	return pow(v_vector[0], 2) + pow(v_vector[1], 2) + pow(v_vector[2], 2);
}

GLfloat glt_get_vector_length(const GLTvector3 v_vector)
{
	return (GLfloat)sqrt(glt_get_vector_length_sqrt(v_vector));
}

void glt_normalize_vector(GLTvector3 v_normal)
{
	GLfloat f_length = 1.0f / glt_get_vector_length(v_normal);
	glt_scale_vector(v_normal, f_length);
}

void glt_draw_torus(GLfloat major_radius, GLfloat minor_radius, GLint num_major, GLint num_minor)
{
	GLTvector3 v_normal;
	double major_step = 2.0f * GLT_PI / num_major;
	double minor_step = 2.0f * GLT_PI / num_minor;
	int i, j;

	//glColor3f(0.0f, 0.0f, 0.5f);

	for(i = 0; i < num_major; ++i)
	{
		double a0 = i * major_step;
		double a1 = a0 + major_step;

		GLfloat x0 = (GLfloat)cos(a0);
		GLfloat y0 = (GLfloat)sin(a0);
		GLfloat x1 = (GLfloat)cos(a1);
		GLfloat y1 = (GLfloat)sin(a1);

		glBegin(GL_TRIANGLE_STRIP);

		for(j = 0; j <= num_minor; ++j)
		{
			double b = j * minor_step;
			GLfloat c = (GLfloat)cos(b);
			GLfloat r = minor_radius * c + major_radius;
			GLfloat z = minor_radius * (GLfloat)sin(b);
			glTexCoord2f((float)(i) / (float)(num_major), (float)(j) / (float)(num_minor));

			v_normal[0] = x0 * c;
			v_normal[1] = y0 * c;
			v_normal[2] = z / minor_radius;

			glVertex3f(x0 * r, y0 * r, z);
			glTexCoord2f((float)(i + 1) / (float)(num_major), (float)(j) / (float)(num_minor));

			v_normal[0] = x1 * c;
			v_normal[1] = y1 * c;
			v_normal[2] = z / minor_radius;

			glNormal3fv(v_normal);
			glVertex3f(x1 * r, y1 * r, z);
		}

		glEnd();
	}
}

void glt_load_identity_matrix(GLTmatrix m)
{
	static GLTmatrix identity = {	1.0f, 0.0f, 0.0f, 0.0f,
					0.0f, 1.0f, 0.0f, 0.0f,
					0.0f, 0.0f, 1.0f, 0.0f,
					0.0f, 0.0f, 0.0f, 1.0f };
	memcpy(m, identity, sizeof(GLTmatrix));
}

void glt_rotation_matrix(float angle, float x, float y, float z, GLTmatrix m_matrix)
{
	float vec_length, sin_save, cos_save, one_minus_cos;
	float xx, yy, zz, xy, yz, zx, xs, ys, zs;

	if(x == 0.0f && y == 0.0f && z == 0.0f)
	{
		glt_load_identity_matrix(m_matrix);
		return;
	}

	vec_length = (float)sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

	x /= vec_length;
	y /= vec_length;
	z /= vec_length;

	sin_save = (float)sin(angle);
	cos_save = (float)cos(angle);
	one_minus_cos = 1.0f - cos_save;

	xx = pow(x, 2);
	yy = pow(y, 2);
	zz = pow(z, 2);
	xy = x * y;
	yz = y * z;
	zx = z * x;
	xs = x * sin_save;
	ys = y * sin_save;
	zs = z * sin_save;

	m_matrix[0] = (one_minus_cos * xx) + cos_save;
	m_matrix[4] = (one_minus_cos * xy) - zs;
	m_matrix[8] = (one_minus_cos * zx) + ys;
	m_matrix[12] = 0.0f;

	m_matrix[1] = (one_minus_cos * xy) + zs;
	m_matrix[5] = (one_minus_cos * yy) + cos_save;
	m_matrix[9] = (one_minus_cos * yz) - xs;
	m_matrix[13] = 0.0f;

	m_matrix[2] = (one_minus_cos * zx) - ys;
	m_matrix[6] = (one_minus_cos * yz) + xs;
	m_matrix[10] = (one_minus_cos * zz) + cos_save;
	m_matrix[14] = 0.0f;

	m_matrix[3] = 0.0f;
	m_matrix[7] = 0.0f;
	m_matrix[11] = 0.0f;
	m_matrix[15] = 1.0f;
}

void draw_torus(void)
{
	GLTmatrix rot_mat, translation_mat, transform_mat;

	static GLfloat y_rot = 0.0f;
	y_rot += 0.5f;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glt_rotation_matrix(DEG2RAD(y_rot), 0.0f, 1.0f, 0.f, transform_mat);
	transform_mat[12] = 0.0f;
	transform_mat[13] = 0.0f;
	transform_mat[14] = -2.5f;

	glLoadMatrixf(transform_mat);

	glt_draw_torus(0.35f, 0.15f, 40, 20);

	glutSwapBuffers();
}

void torus_reshape(int w, int h)
{
	GLfloat f_aspect;

	if(h == 0)
		h = 1;

	glViewport(0, 0, w, h);

	f_aspect = (GLfloat)w / (GLfloat)h;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(35.0f, f_aspect, 1.0f, 50.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void torus_timer(int val)
{
	glutPostRedisplay();
	glutTimerFunc(33, torus_timer, 1);
}
