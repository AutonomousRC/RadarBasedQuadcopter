#include "bluestein_fft.h"
#include "complex_math.h"
#include <stdbool.h>
#include <math.h>

void bluestein_first(complex *x, int init, double *costab, double *sintab, complex *y)
{
	int i, j, ix, ju, iy, iheight, istart, ihi;
	double temp_re, temp_im, twid_re, twid_im;
	bool tst;

	for(i = 0; i < 2048; i++)
	{
		y[i].re = 0.0;
		y[i].im = 0.0;
	}

	ix = init;
	ju = 0;
	iy = 0;

	for(i = 0; i < 799; i++)
	{
		y[iy] = x[ix];
		iy = 2048;
		tst = true;
		while (tst)
		{
			iy >>= 1;
			ju ^= iy;
			tst = ((ju & iy) == 0);
		}

		iy = ju;
		ix++;
	}

	y[iy] = x[ix];

	for(i = 0; i <= 2047; i += 2)
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
	ju = 512;
	iheight = 2045;

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

		for(j = ju; j < 1024; j += ju)
		{
			twid_re = costab[j];
			twid_im = sintab[j];
			i = istart;
			ihi = istart + iheight;

			while(i < ihi)
			{
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
}

void bluestein_second(complex *x, double *costab, double *sintab, complex *y)
{
}

void set_bluestein(complex *wwc)
{
	int i, rt, idx, y;
	double nt_im;

	idx = 798;
	rt = 0;

	wwc[779].re = 1.0;
	wwc[779].im = 0.0;

	for(i = 0; i < 799; i++)
	{
		y = ((i + 1) << 1) - 1;

		if(1600 - rt <= y)
			rt = (y + rt) - 1600;
		else
			rt += y;

		nt_im = -M_PI * (double)rt / 800.0;
		wwc[idx].re = cos(nt_im);
		wwc[idx].im = -sin(nt_im);
		idx--;
	}

	idx = 0;

	for(i = 798; i >= 0; i--)
	{
		wwc[i + 800] = wwc[idx];
		idx++;
	}
}

void fft(complex *x, complex *y)
{
	complex wwc[1599] = {0};
	complex res[800] = {0};

	complex fy[2048] = {0};
	complex fv[2048] = {0};

	int i, xidx;
	double t, step, fy_re;

	double costab[1025] = {0};
	double sintab[1025] = {0};
	double sintab_inv[1025] = {0};

	step = 2 * M_PI / 2048;

	for(i = 0; i < 1025; i++)
	{
		costab[i] = cos(t);
		sintab[i] = sin(t);
		sintab_inv[i] = -sintab[i];
		t += step;
	}

	set_bluestein(wwc);

	xidx = 0;

	for(i = 0; i < 800; i++)
	{
		res[i].re = wwc[i + 799].re * x[xidx].re + wwc[i + 799].im * x[xidx].im;
		res[i].im = wwc[i + 799].re * x[xidx].im - wwc[i + 799].im * x[xidx].re;
		xidx++;
	}

	bluestein_first(res, 0, costab, sintab, fy);
	bluestein_second(fy, costab, sintab, fv);
}
