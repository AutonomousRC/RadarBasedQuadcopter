#include "radar_basic.h"
#include "n2_fft.h"
#include "ogl_helper.h"
#include "complex_math.h"
#include "math_tech.h"
#include "init_require_data.h"
#include "hanning_window.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define	SAMPLE		800

void set_hanning(double *data, int num)
{
	int i;

	for(i = 0; i < num; i++)
		data[i] = 0.5 * (1 - cos(2 * M_PI * i / (num - 1)));
}

void draw_hanning_windowed_fft(void)
{
	int i, cache = 0;
	double max, min, cx, cy, epsilon = 0.001;
	comp freq[256] = {0};
	double sample[32] = {0};
	double freq_amp[256] = {0};
	double tmp[256] = {0};

	set_hanning(sample, 32);

	n2_fft(sample, freq, 32, 256);
	complex_abs(freq, freq_amp, 256);

	memcpy(&tmp[0], &freq_amp[0], sizeof(double) * 256);
	max = find_max(tmp, 256);

	divide_vec(freq_amp, max, 256);

	volt_base_log(freq_amp, epsilon, 256);
	memcpy(&tmp[0], &freq_amp[0], sizeof(double) * 256);
	min = find_min(tmp, 256);

	glBegin(GL_LINES);

        for(i = 0; i < 256; i++)
        {
                if(cache)
                {
                        glVertex2f(-130 + i * (260.0 / 256.0), 80.0 + freq_amp[i] * (-160.0 / min));
                        glVertex2f(-130 + cx * (260.0 / 256.0), 80.0 + cy * (-160.0 / min));
                }

                cache = 1;
                cx = i;
                cy = freq_amp[i];
        }

        glEnd();
}

void hanning_window(void)
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

	draw_hanning_windowed_fft();

	glutSwapBuffers();
}
