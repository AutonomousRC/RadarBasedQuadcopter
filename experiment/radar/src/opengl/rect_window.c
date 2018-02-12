#include "radar_basic.h"
#include "n2_fft.h"
#include "ogl_helper.h"
#include "complex_math.h"
#include "math_tech.h"
#include "init_require_data.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define	SAMPLE		800

#if 0
// first case of sample num = 800
void draw_lfm_fft(void)
{
	int i, cache = 0;
	double sig[SAMPLE] = {0};
	double freq[SAMPLE] = {0};
	double cx, cy;

	init_require_data(sig, SAMPLE);

	//p_arr(sig, 800);

#if 1
	linear_slice(-50000000, 50000000, SAMPLE, freq);

	//p_arr(freq, 800);

	glBegin(GL_LINES);

	for(i = 0; i < SAMPLE; i++)
	{
		if(cache)
		{
			glVertex2f(freq[i] * (130.0 / 50000000.0), -80.0 + sig[i] * 3.0);
			glVertex2f(cx * (130.0 / 50000000.0), -80.0 + cy * 3.0);
		}

		cache = 1;
		cx = freq[i];
		cy = sig[i];
	}

	glEnd();
#endif
}
#endif

void draw_rect_windowed_fft(void)
{
	int i, cache = 0;
	double max, min, cx, cy, epsilon = 0.001;
	complex freq[256] = {0};
	double sample[32] = {0};
	double freq_amp[256] = {0};
	double tmp[256] = {0};

	for(i = 0; i < 32; i++)
		sample[i] = 1.0;

	n2_fft(sample, freq, 32, 256);
	complex_abs(freq, freq_amp, 256);

#if 0
	for(i = 0; i < 256; i++)
		printf("freq_amp[%d] = %lf\n", i, freq_amp[i]);
#endif

	memcpy(&tmp[0], &freq_amp[0], sizeof(double) * 256);
	max = find_max(tmp, 256);

#if 0
	for(i = 0; i < 256; i++)
		printf("freq_amp[%d] = %lf\n", i, freq_amp[i]);
#endif

#if 0
	printf("max = %lf\n", max);
#endif

	divide_vec(freq_amp, max, 256);

#if 0
	for(i = 0; i < 256; i++)
		printf("freq_amp[%d] = %lf\n", i, freq_amp[i]);
#endif

	volt_base_log(freq_amp, epsilon, 256);
	memcpy(&tmp[0], &freq_amp[0], sizeof(double) * 256);
	min = find_min(tmp, 256);
	//max = find_max(tmp, 256);
	//printf("min = %lf\n", min);
	//printf("max = %lf\n", max);

#if 0
	for(i = 0; i < 256; i++)
		printf("freq_amp[%d] = %lf\n", i, freq_amp[i]);
#endif

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

void rect_window(void)
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

	//draw_lfm_fft();
	draw_rect_windowed_fft();

	glutSwapBuffers();
}
