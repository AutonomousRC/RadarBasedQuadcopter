#include "n2_fft.h"
#include "complex_math.h"
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#if 0
int calc_align_idx(int num)
{
	int i, tmp = 0;

	for(i = 1; ; i++)
	{
		tmp = 1 << i;
		if(tmp > num)
			return i;
	}
}
#endif

// first case x_num 32, y_num 256
void n2_fft(double *x, complex *y, int x_num, int y_num)
{
	bool tst;
	complex res[256] = {0};
	double temp_re, temp_im, twid_re, twid_im;
	int i, ix, iy, j, ju, istart, iheight, ihi;

	double step, t = 0;

	double costab[129] = {0};
	double sintab[129] = {0};

	int yn = y_num;
	int ynm1 = yn - 1;
	int ynm3 = yn - 3;
	int ynd2 = yn / 2;
	int ynd2p1 = ynd2 + 1;
	int ynd4 = yn / 4;

	int xn = x_num;
	int xnm1 = xn - 1;

	//step = 2 * M_PI / 256;
	step = 2 * M_PI / yn;

	//for(i = 0; i < 129; i++)
	for(i = 0; i < ynd2p1; i++)
	{
		costab[i] = cos(t);
		sintab[i] = -sin(t);
		t += step;

#if 0
		printf("costab[%d] = %lf\n", i, costab[i]);
		printf("sintab[%d] = %lf\n", i, sintab[i]);
#endif
	}

	ix = iy = ju = 0;

	//for(i = 0; i < 31; i++)
	for(i = 0; i < xnm1; i++)
	{
		res[iy].re = x[ix];
		res[iy].im = 0.0;
		//iy = 256;
		iy = yn;
		tst = true;

		while(tst)
		{
			iy >>= 1;
			ju ^= iy;
			tst = ((ju & iy) == 0);
		}

		iy = ju;
		ix++;
	}

	res[iy].re = x[ix];
	res[iy].im = 0.0;

	for(i = 0; i <= ynm1; i += 2)
	{
		temp_re = res[i + 1].re;
		temp_im = res[i + 1].im;
		res[i + 1].re = res[i].re - res[i + 1].re;
		res[i + 1].im = res[i].im - res[i + 1].im;
		res[i].re += temp_re;
		res[i].im += temp_im;
	}

	iy = 2;
	ix = 4;
	//ju = 64;
	ju = ynd4;
	//iheight = 253;
	iheight = ynm3;

	while(ju > 0)
	{
		for(i = 0; i < iheight; i += ix)
		{
			temp_re = res[i + iy].re;
			temp_im = res[i + iy].im;
			res[i + iy].re = res[i].re - temp_re;
			res[i + iy].im = res[i].im - temp_im;
			res[i].re += temp_re;
			res[i].im += temp_im;
		}

		istart = 1;

		//for(j = ju; j < 128; j += ju)
		for(j = ju; j < ynd2; j += ju)
		{
			twid_re = costab[j];
			twid_im = sintab[j];

			i = istart;
			ihi = istart + iheight;

			while(i < ihi)
			{
				temp_re = twid_re * res[i + iy].re - twid_im * res[i + iy].im;
				temp_im = twid_re * res[i + iy].im + twid_im * res[i + iy].re;
				res[i + iy].re = res[i].re - temp_re;
				res[i + iy].im = res[i].im - temp_im;
				res[i].re += temp_re;
				res[i].im += temp_im;
				i += ix;
			}

			istart++;
		}

		ju /= 2;
		iy = ix;
		ix <<= 1;
		iheight -= iy;
	}

	memcpy(&y[0], &res[0], sizeof(complex) << 8);
}

#if 0
void fft_shift(double *x)
{
	int i, i1, i2, ib, a, k, dim, vlend2, npages;
	double xtmp;

	for(dim = 0; dim < 2; dim++)
	{
		a = 1 + 799 * dim;

		if (a <= 1)
			;
		else
		{
			vlend2 = a / 2;
			npages = 1;
			k = dim + 2;

			while(k < 3)
			{
				npages *= 800;
				k = 3;
			}

			if(vlend2 << 1 == a)
			{
				i2 = 0;

				for(i = 1; i <= npages; i++)
				{
					i1 = i2;
					i2 += a;
					ib = i1 + vlend2;

					for(k = 1; k <= vlend2; k++)
					{
						xtmp = x[i1];
						x[i1] = x[ib];
						x[ib] = xtmp;
						i1++;
						ib++;
					}
				}
			}
			else
			{
				i2 = 0;

				for(i = 1; i <= npages; i++)
				{
					i1 = i2;
					i2 += a;
					ib = i1 + vlend2;
					xtmp = x[ib];

					for(k = 1; k <= vlend2; k++)
					{
						x[ib] = x[i1];
						x[i1] = x[ib + 1];
						i1++;
						ib++;
					}

					x[ib] = xtmp;
				}
			}
		}
	}
}
#endif
