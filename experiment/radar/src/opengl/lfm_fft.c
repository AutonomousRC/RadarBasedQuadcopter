#include "radar_basic.h"
#include "lfm_fft.h"
#include "ogl_helper.h"
#include "math_tech.h"
#include "init_require_data.h"
#include <stdlib.h>
#include <math.h>

void draw_lfm_fft(void)
{
	int i;
	double sig[800] = {0};
	double freq[800] = {0};

	init_require_data(sig);

#if 1
	linear_slice(-50000000, 50000000, 800, freq);

	//p_arr(freq, 800);

	glBegin(GL_LINES);

	for(i = 0; i < 800; i++)
	{
		glVertex2f(freq[i] * (130.0 / 50000000.0), sig[i] * (80.0 / 1000.0));
	}

	glEnd();
#endif
}

void lfm_fft(void)
{
	static char label[100];
	float vert[11] = {0};
	float hori[5] = {0};
	int i;

	init_screen();
	rect_screen();

	slc(-130, 130, 8, vert);
	slc(-80, 80, 6, hori);

	for(i = 0; i < 8; i++)
	{
		glBegin(GL_LINE_LOOP);
		glVertex3f(vert[i], 80.0, 0.0);
		glVertex3f(vert[i], -80.0, 0.0);
		glEnd();
	}

	for(i = 0; i < 6; i++)
	{
		glBegin(GL_LINE_LOOP);
		glVertex3f(130.0, hori[i], 0.0);
		glVertex3f(-130.0, hori[i], 0.0);
		glEnd();
	}

	glColor3f(1, 0, 0);

	draw_lfm_fft();

	glutSwapBuffers();
}
