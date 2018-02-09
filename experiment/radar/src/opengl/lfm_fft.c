#include "radar_basic.h"
#include "lfm_fft.h"
#include "ogl_helper.h"
#include "math_tech.h"
#include <stdlib.h>
#include <math.h>

#define SLICE                   (1024)
#define HALF_SLICE              (SLICE >> 1)
#define QUAD_SLICE              (SLICE >> 2)
#define CALC_ORDER              ((HALF_SLICE) + 1)
#define CALC_HEIGHT             (SLICE - 3)

#if 0

#define	SLICE		(800)
#define HALF_SLICE	(SLICE >> 1)
#define	CALC_ORDER	((HALF_SLICE) + 1)

void find_frequency(c X[SLICE], double R[CALC_ORDER], double f[HALF_SLICE])
{
        double P2[SLICE];
        c y[SLICE];
        int k;
        for (k = 0; k < SLICE; k++) {
#if 0
                if (X[k].im == 0.0) {
                        y[k].re = X[k].re / 256.0;
                        y[k].im = 0.0;
                } else if (X[k].re == 0.0) {
                        y[k].re = 0.0;
                        y[k].im = X[k].im / 256.0;
                } else {
#endif
                y[k].re = X[k].re / SLICE;
                y[k].im = X[k].im / SLICE;
//              }

                P2[k] = rt_hypotd_snf(y[k].re, y[k].im);
                //printf("P2[%d] = %lf\n", k, P2[k]);
        }

        memcpy(&R[0], &P2[0], CALC_ORDER * sizeof(double));
        for (k = 0; k < HALF_SLICE - 1; k++) {
                R[1 + k] = 2.0 * P2[1 + k];
                //printf("R[%d] = %lf\n", k + 1, R[k + 1]);
        }

        for(k = 0; k < CALC_ORDER; k++)
        {
                f[k] = SAMPLE_FREQ * k / SLICE;
                //printf("f[%d] frequency = %lf\n", k, f[k]);
        }
}

#endif

void draw_lfm_fft(void)
{
	int k, nscat = 2, cache = 0, n;
	double taup = 10 * pow(10, -6);
	double b = 40 * pow(10, 6);
	double rrec = 50;
	double scat_range[3] = {15, 25};
	double scat_rcs[3] = {1, 2};
	double eps = pow(10, -16);
	double time_b_prod = b * taup;

	c y[800] = {0};
	int ix = 0, ju = 0, iy = 0, tst, iheight, istart, ihi, i, j;
	double temp_re, temp_im, twid_re, twid_im;
	double dv0[800] = {0};
	double dv1[800] = {0};
	double rf[800] = {0};
	double f[800] = {0};

	double step, cx, cy;
	double t[801] = {0};
	double signal[801] = {0};

	double sampling_interval;

	if(time_b_prod < 5)
	{
		printf("Time Bandwidth Term is too small\n");
		exit(-1);
	}

	n = 2 * taup * b + 1;

	//printf("n = %d\n", n);

	linear_slice(-taup / 2.0, taup / 2.0, n, t);

	//p_arr(t, 800);

	for(k = 0; k < 800; k++)
		signal[k] = sin(M_PI * (b / taup) * pow(t[k], 2));

	p_arr(signal, 800);

	//slice_section(0.0, tau, ts, t);
	//calc_period(&freq, &period);
	//step = 0.001;

#if 0
	for(i = 0; i < 401; i++)
	{
		dv0[i] = cos(t[i]);
		dv1[i] = -sin(t[i]);
	}

	for(i = 0; i < 799; i++)
	{
		y[iy].re = signal[i];
		y[iy].im = 0.0;

		iy = 800;
		tst = 1;

		while(tst)
		{
			iy >>= 1;
			ju ^= iy;
			tst = ((ju & iy) == 0);
		}

		iy = ju;
		ix++;
	}

	y[iy].re = signal[ix];
	y[iy].im = 0.0;

	for(i = 0; i <= 800 - 1; i += 2)
	{
		temp_re = y[i + 1].re;
		temp_im = y[i + 1].im;
		y[i + 1].re = y[i].re - y[i + 1].re;
		y[i + 1].im = y[i].im - y[i + 1].im;
		y[i].re += temp_re;
		y[i].im += temp_im;
	}

	iy = 2;
	ix = 4;
	ju = 200;
	iheight = 797;

	while(ju > 0)
	{
		for(i = 0; i < iheight; i += ix)
		{
			temp_re = y[i + iy].re;
			temp_im = y[i + iy].im;
			y[i + iy].re = y[i].re - temp_re;
			y[i + iy].im = y[i].im - temp_im;
			y[i].re += temp_re;
			y[i].im += temp_im;
		}

		istart = 1;

		for(j = ju; j < 400; j += ju)
		{
			twid_re = dv0[j];
			twid_im = dv1[j];
			i = istart;
                        ihi = istart + iheight;
                        while (i < ihi) {
                                temp_re = twid_re * y[i + iy].re - twid_im * y[i + iy].im;
                                temp_im = twid_re * y[i + iy].im + twid_im * y[i + iy].re;
                                y[i + iy].re = y[i].re - temp_re;
                                y[i + iy].im = y[i].im - temp_im;
                                y[i].re += temp_re;
                                y[i].im += temp_im;
                                i += ix;
                        }

                        istart++;
		}

		ju /= 2;
                iy = ix;
                ix <<= 1;
                iheight -= iy;
	}

	//if(step > 1)
	//	step = 0.0;

	for(i = 0; i < 401; i++)
	{
		printf("tx = %lf\n", t[i] * 100000000);
		printf("y = %lf\n", y[i].re * 16);

		glBegin(GL_LINE_STRIP);
		glVertex2f(t[i] * 100000000, y[i].re * 16);
		glVertex2f(t[i] * 100000000, 0);
		glEnd();
	}

	printf("sampling interval = %.12lf\n", 1 / 2.5 / b);
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
