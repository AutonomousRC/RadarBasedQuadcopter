#include "lfm_visual.h"
#include "ogl_helper.h"
#include "math_tech.h"
#include <math.h>

void calc_period(float *freq, float *period)
{
	*period = 1 / (*freq);
}

void draw_lfm_sin_signal(void)
{
	double lfm_bandwidth = 15 * pow(10, 6);
	double tau = pow(10, -6);
	double ts = pow(10, -9);
	double beta = lfm_bandwidth / tau;

	double t[1002] = {0};

	double y, step, cx, cy, period, freq = 15;
	int i, cache = 0;

	slice_section(0.0, tau, ts, t);


	calc_period(&freq, &period);

	step = 0.001;

	if(step > 1)
		step = 0.0;

	glBegin(GL_LINES);
	
	for(i = 0; i < 1001; i++)
	{
		y = sin(beta * M_PI * pow(t[i], 2));
		//printf("y = %lf\n", y * 80);
		printf("x = %lf\n", t[i] * 260000000);

		if(cache)
		{
			glVertex2f(-130 + cx * 260000000, cy * 80);
			glVertex2f(-130 + t[i] * 260000000, y * 80);
		}

		cache = 1;
		cx = t[i];
		cy = y;
	}

	glEnd();
}

void lfm_visual(void)
{
	static char label[100];
	float vert[11] = {0};
	float hori[5] = {0};
	int i;

	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glColor3f(0, 1, 0);

	glBegin(GL_LINE_LOOP);
	glVertex3f(130.0, -80.0, 0.0);
	glVertex3f(-130.0, -80.0, 0.0);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(130.0, 80.0, 0.0);
	glVertex3f(-130.0, 80.0, 0.0);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(-130.0, 80.0, 0.0);
	glVertex3f(-130.0, -80.0, 0.0);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(130.0, 80.0, 0.0);
	glVertex3f(130.0, -80.0, 0.0);
	glEnd();

#if 0
	glBegin(GL_LINE_LOOP);
	glVertex3f(0.0, 80.0, 0.0);
	glVertex3f(0.0, -80.0, 0.0);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(130.0, 0.0, 0.0);
	glVertex3f(-130.0, 0.0, 0.0);
	glEnd();
#endif

	slc(-130, 130, 10, vert);
	slc(-80, 80, 4, hori);

	for(i = 0; i < 10; i++)
	{
		glBegin(GL_LINE_LOOP);
		glVertex3f(vert[i], 80.0, 0.0);
		glVertex3f(vert[i], -80.0, 0.0);
		glEnd();
	}

	for(i = 0; i < 4; i++)
	{
		glBegin(GL_LINE_LOOP);
		glVertex3f(130.0, hori[i], 0.0);
		glVertex3f(-130.0, hori[i], 0.0);
		glEnd();
	}

	glColor3f(1, 0, 0);

	draw_lfm_sin_signal();

	glutSwapBuffers();
}
