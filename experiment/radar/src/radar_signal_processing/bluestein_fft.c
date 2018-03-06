#include "bluestein_fft.h"
#include "complex_math.h"
//
#include "init_require_data.h"
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

void bluestein_first(comp *x, int init, double *costab, double *sintab, comp *y, int data_num)
{
	int i, j, ix, ju, iy, iheight, istart, ihi;
	double temp_re, temp_im, twid_re, twid_im;
	bool tst;

#if 1
	int dn = data_num;
	int cdn = 1 << calc_align_idx(dn);
	int cdn2 = cdn * 2;
	int dnm1 = dn - 1;
	int cdn2m1 = cdn2 - 1;
	int cdn2m3 = cdn2 - 3;
	int cdnd2 = cdn / 2;
#endif

	//for(i = 0; i < 2048; i++)
	for(i = 0; i < cdn2; i++)
	{
		y[i].re = 0.0;
		y[i].im = 0.0;
	}

	ix = init;
	ju = 0;
	iy = 0;

	//for(i = 0; i < 799; i++)
	for(i = 0; i < dnm1; i++)
	{
		y[iy] = x[ix];
		//iy = 2048;
		iy = cdn2;
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

	//for(i = 0; i <= 2047; i += 2)
	for(i = 0; i <= cdn2m1; i += 2)
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
	/*
	ju = 512;
	iheight = 2045;
	*/
	ju = cdnd2;
	iheight = cdn2m3;

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

		//for(j = ju; j < 1024; j += ju)
		for(j = ju; j < cdn; j += ju)
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

void bluestein_second(comp *x, double *costab, double *sintab, comp *y, int data_num)
{
	bool tst;
	double temp_re, temp_im, twid_re, twid_im;
	int i, j, ix, iy, ju, iheight, istart, ihi;

#if 1
	int dn = data_num;			// first case 1599
	int cdn = 1 << calc_align_idx(dn);	// first case 2048
	int dnm1 = dn - 1;			// first case 1598
	int cdnm1 = cdn - 1;			// first case 2047
	int cdnm3 = cdn - 3;			// first case 2045
	int cdnd2 = cdn / 2;			// first case 1024
	int cdnd4 = cdn / 4;			// first case 512

	//printf("dn = %d, cdn = %d, dnm1 = %d, cdnm1 = %d, cdnm3 = %d, cdnd2 = %d, cdnd4 = %d\n", dn, cdn, dnm1, cdnm1, cdnm3, cdnd2, cdnd4);
#endif

	//for(i = 0; i < 2048; i++)
	for(i = 0; i < cdn; i++)
	{
		y[i].re = 0.0;
		y[i].im = 0.0;
	}

	ix = 0;
	ju = 0;
	iy = 0;

	//for(i = 0; i < 1598; i++)
	for(i = 0; i < dnm1; i++)
	{
		y[iy] = x[ix];
		//iy = 2048;
		iy = cdn;
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

	//for (i = 0; i <= 2047; i += 2)
	for (i = 0; i <= cdnm1; i += 2)
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
	//ju = 512;
	ju = cdnd4;
	//iheight = 2045;
	iheight = cdnm3;

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

		//for(j = ju; j < 1024; j += ju)
		for(j = ju; j < cdnd2; j += ju)
		{
			twid_re = costab[j];
			twid_im = sintab[j];
			i = istart;
			ihi = istart + iheight;

			while (i < ihi)
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

void bluestein_third(comp *x, double *costab, double *sintab, comp *y, int data_num)
{
	int i, ix, iy, j, ju, iheight, istart, ihi;
	double temp_re, temp_im, twid_re, twid_im;
	bool tst;

#if 1
	int dn = data_num;	// first case 2048
	int dnm1 = dn - 1;	// first case 2047
	int dnd2 = dn / 2;	// first case 1024
	int dnd4 = dn / 4;	// first case 512
	int dnm3 = dn - 3;	// first case 2045

	//printf("dn = %d, dnm1 = %d\n", dn, dnm1);
#endif

	ix = 0;
	ju = 0;
	iy = 0;

	//for(i = 0; i < 2047; i++)
	for(i = 0; i < dnm1; i++)
	{
		y[iy] = x[ix];
		//iy = 2048;
		iy = dn;
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

	y[iy] = x[ix];

	//for(i = 0; i <= 2047; i += 2)
	for(i = 0; i <= dnm1; i += 2)
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
	//ju = 512;
	ju = dnd4;
	//iheight = 2045;
	iheight = dnm3;

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

		//for(j = ju; j < 1024; j += ju)
		for(j = ju; j < dnd2; j += ju)
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

	//for(iy = 0; iy < 2048; iy++)
	for(iy = 0; iy < dn; iy++)
	{
		y[iy].re *= sin(2 * M_PI / dn) / (2 * M_PI);
		y[iy].im *= sin(2 * M_PI / dn) / (2 * M_PI);
	/*
		y[iy].re *= sin(2 * M_PI / 2048) / (2 * M_PI);
		y[iy].im *= sin(2 * M_PI / 2048) / (2 * M_PI);
	*/
		//y[iy].re *= 0.00048828125;
		//y[iy].im *= 0.00048828125;
	}
}

void set_bluestein(comp *wwc, int data_num)
{
	int i, rt, idx, y;
	double nt_im;

#if 1
	int dn = data_num;
	int dn2 = dn * 2;
	int dnm1 = dn - 1;
	int dnm2 = dn - 2;
#endif

	//idx = 798;
	idx = dnm2;
	rt = 0;

	/*
	wwc[799].re = 1.0;
	wwc[799].im = 0.0;
	*/
	wwc[dnm1].re = 1.0;
	wwc[dnm1].im = 0.0;

	for(i = 0; i < dnm1; i++)
	{
		y = ((i + 1) << 1) - 1;

		//if(1600 - rt <= y)
		if(dn2 - rt <= y)
			//rt = (y + rt) - 1600;
			rt = (y + rt) - dn2;
		else
			rt += y;

		//nt_im = -M_PI * (double)rt / 800.0;
		nt_im = -M_PI * (double)rt / (double)dn;
		wwc[idx].re = cos(nt_im);
		wwc[idx].im = -sin(nt_im);

#if 0
		printf("wwc[%d].re = %lf\n", idx, wwc[idx].re);
		printf("wwc[%d].im = %lf\n", idx, wwc[idx].im);
#endif

		idx--;
	}

	//printf("wwc[799].re = %lf\n", wwc[799].re);

	idx = 0;

	//for(i = 798; i >= 0; i += -1)
	for(i = dnm2; i >= 0; i += -1)
	{
		//wwc[i + 800] = wwc[idx];
		wwc[i + dn] = wwc[idx];
		idx++;
	}
}

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

void bluestein_fft(comp *x, comp *y, int y_data_num)
{
	/* Can't Use this method because of heap corruption problem */

	/* dn2m1 = wwc = y_data_num * 2 - 1
	   res = y_data_num
	   cdn2 = fy = calc_align_idx(y_data_num) * 2
	         fv = calc_align_idx(y_data_num) * 2
	   cdn = costab = calc_align_idx(y_data_num)
	         sintab = calc_align_idx(y_data_num)
	         sintab_inv = calc_align_idx(y_data_num)
	   step = 2 * M_PI / (calc_align_idx(y_data_num) * 2)
	   tmp = step * calc_align_idx(y_data_num)
	   dnm1 = y_data_num - 1
	   cdnp1 = for loop = calc_align_idx(y_data_num) + 1 */

#if 1
	int dn = y_data_num;
	int dn2m1 = dn * 2 - 1;
	int cdn = 1 << calc_align_idx(dn);
	int cdn2 = cdn * 2;
	int dnm1 = dn - 1;
	int cdnp1 = cdn + 1;
#endif
#if 0
	comp *wwc = NULL;
	comp *res = NULL;
	comp *fy = NULL;
	comp *fv = NULL;
#endif
	comp wwc[1599] = {0};
	comp res[800] = {0};

	comp fy[2048] = {0};
	comp fv[2048] = {0};

	int i, xidx;
	double t = 0, step, fy_re;
	double tmp;

	double costab[1025] = {0};
	double sintab[1025] = {0};
	double sintab_inv[1025] = {0};

	/* Heap Corruption */
#if 0
	complex_arr_alloc(&wwc, dn2m1);
	complex_arr_alloc(&res, dn);
	complex_arr_alloc(&fy, cdn2);
	complex_arr_alloc(&fv, cdn2);
#endif

	//step = 2 * M_PI / 2048;
	step = 2 * M_PI / cdn2;
	//tmp = step * 1024;
	tmp = step * cdn;

	//printf("dn = %d, dn2m1 = %d, cdn = %d, cdn2 = %d, dnm1 = %d, cdnp1 = %d\n", dn, dn2m1, cdn, cdn2, dnm1, cdnp1);

	// for(i = 0; i < 1025; i++)
	for(i = 0; i < cdnp1; i++)
	{
		costab[i] = cos(t);
		sintab[i] = -sin(t);
		sintab_inv[i] = -sin(tmp + t);
		t += step;
#if 0
		printf("costab[%d] = %lf\n", i, costab[i]);
		printf("sintab[%d] = %lf\n", i, sintab[i]);
		printf("sintab_inv[%d] = %lf\n", i, sintab_inv[i]);
#endif
	}

	set_bluestein(wwc, dn);

	//print_complex(wwc, 1599);

	/* y_data_num
	*/

	xidx = 0;

	//for(i = 0; i < 800; i++)
	for(i = 0; i < dn; i++)
	{
		/*
		res[i].re = wwc[i + 799].re * x[xidx].re + wwc[i + 799].im * x[xidx].im;
		res[i].im = wwc[i + 799].re * x[xidx].im - wwc[i + 799].im * x[xidx].re;
		*/
		res[i].re = wwc[i + dnm1].re * x[xidx].re + wwc[i + dnm1].im * x[xidx].im;
		res[i].im = wwc[i + dnm1].re * x[xidx].im - wwc[i + dnm1].im * x[xidx].re;
		xidx++;
	}

	//print_complex(res, 800);

#if 0
	printf("wwc[799].re = %lf\n", wwc[799].re);
	printf("x[0].re = %lf\n", x[0].re);
	printf("wwc[799].im = %lf\n", wwc[799].im);
	printf("x[0].im = %lf\n", x[0].im);
#endif

	bluestein_first(res, 0, costab, sintab, fy, dn);

	//print_complex(fy, 2048);

	bluestein_second(wwc, costab, sintab, fv, dn2m1);

	//print_complex(fv, 2048);

	//for(xidx = 0; xidx < 2048; xidx++)
	for(xidx = 0; xidx < cdn2; xidx++)
	{
		fy_re = fy[xidx].re; 
		fy[xidx].re = fy[xidx].re * fv[xidx].re - fy[xidx].im * fv[xidx].im;
		fy[xidx].im = fy_re * fv[xidx].im + fy[xidx].im * fv[xidx].re;
	}

	//print_complex(fy, 2048);

	bluestein_third(fy, costab, sintab_inv, fv, cdn2);

	//print_complex(fv, 2048);

	xidx = 0;

	//for(i = 0; i < 800; i++)
	for(i = 0; i < dn; i++)
	{
		/*
		res[xidx].re = wwc[i + 799].re * fv[i + 799].re + wwc[i + 799].im * fv[i + 799].im;
		res[xidx].im = wwc[i + 799].re * fv[i + 799].im - wwc[i + 799].im * fv[i + 799].re;
		*/
		res[xidx].re = wwc[i + dnm1].re * fv[i + dnm1].re + wwc[i + dnm1].im * fv[i + dnm1].im;
		res[xidx].im = wwc[i + dnm1].re * fv[i + dnm1].im - wwc[i + dnm1].im * fv[i + dnm1].re;
		xidx++;
	}

	//print_complex(res, 800);

	//memcpy(&y[0], &res[0], 800U * sizeof(comp));
	memcpy(&y[0], &res[0], dn * sizeof(comp));
}

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
