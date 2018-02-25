#include "zeroth_modify_bessel.h"
#include "math_tech.h"
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

// bool = bool
// MIN_int32_T = -2147483648

double rtNaN;
double rtInf;
double rtMinusInf;

float rtNaNF;
float rtInfF;
float rtMinusInfF;

void zeroth_modify_bessel(double *data, complex *z, int num)
{
	double y[32];
	complex zd, b_w;
	int k, b_k, ierr, unused;

	memcpy(&y[0], &data[0], sizeof(double) << 5);

#if 0
	for(k = 0; k < 32; k++)
		printf("data[%d] = %lf\n", k, data[k]);
#endif

	for (k = 1; k < 33; k++)
	{
		b_k = k;
		zd.re = y[b_k - 1];
		zd.im = 0.0;

		cbesi(zd, 0.0, 1, &b_w, &unused, &ierr);

#if 1
		printf("%d. b_w.re = %lf\n", k, b_w.re);
		printf("%d. b_w.im = %lf\n", k, b_w.im);
#endif

		if (y[b_k - 1] > 0.0)
			b_w.im = 0.0;

		z[b_k - 1] = b_w;
	}
}

int casyi(complex z, double fnu, int kode, int nin, complex *y, double rl, double tol, double elim)
{
	int i;
	int nz;
	int n;
	double x;
	complex ak1;
	double brm;
	double acz;
	double s;
	double cz_re;
	double cz_im;
	double dnu2;
	double fdn;
	double ez_re;
	double ez_im;
	double aez;
	double b_dnu2;
	int jl;
	double p1_re;
	int inu;
	double bk;
	bool exitg1;
	double sqk;
	double b_atol;
	double sgn;
	double cs1_re;
	double cs1_im;
	double cs2_re;
	double cs2_im;
	double ak;
	double aa;
	double bb;
	double dk_re;
	double dk_im;
	bool errflag;
	bool exitg2;
	double b_cz_re;
	complex test = {2, 2};

	printf("casyi\n");

	if(nin <= 1)
		n = nin;
	else
		n = 1;

	nz = 0;
	x = z.re;
	if(z.im == 0.0)
	{
		ak1.re = 0.15915494309189535 / z.re;
		ak1.im = 0.0;
	}
	else if(z.re == 0.0)
	{
		ak1.re = 0.0;
		ak1.im = -(0.15915494309189535 / z.im);
	}
	else
	{
		brm = fabs(z.re);
		acz = fabs(z.im);
		if(brm > acz)
		{
			s = z.im / z.re;
			acz = z.re + s * z.im;
			ak1.re = (0.15915494309189535 + s * 0.0) / acz;
			ak1.im = (0.0 - s * 0.15915494309189535) / acz;
		}
		else if(acz == brm)
		{
			if(z.re > 0.0)
				acz = 0.5;
			else
				acz = -0.5;

			if (z.im > 0.0)
				dnu2 = 0.5;
			else 
				dnu2 = -0.5;

			ak1.re = (0.15915494309189535 * acz + 0.0 * dnu2) / brm;
			ak1.im = (0.0 * acz - 0.15915494309189535 * dnu2) / brm;
		}
		else
		{
			s = z.re / z.im;
			acz = z.im + s * z.re;
			ak1.re = s * 0.15915494309189535 / acz;
			ak1.im = (s * 0.0 - 0.15915494309189535) / acz;
		}
	}

	b_sqrt(&ak1);
	b_sqrt(&test);
	printf("test.re = %lf\n", test.re);
	printf("test.im = %lf\n", test.im);

	if(kode == 2)
	{
		cz_re = 0.0;
		cz_im = z.im;
		acz = 0.0;
	}
	else
	{
		cz_re = z.re;
		cz_im = z.im;
		acz = z.re;
	}

	if(fabs(acz) > elim)
	{
		nz = -1;
		y->re = rtNaN;
		y->im = 0.0;
	}
	else
	{
		dnu2 = fnu + fnu;
		if(rtIsInf(cz_im) && rtIsInf(cz_re) && (cz_re < 0.0))
		{
			cz_re = 0.0;
			cz_im = 0.0;
		}
		else
		{
			acz = exp(cz_re / 2.0);
			cz_re = acz * (acz * cos(cz_im));
			cz_im = acz * (acz * sin(cz_im));
		}

		acz = ak1.re;
		ak1.re = ak1.re * cz_re - ak1.im * cz_im;
		ak1.im = acz * cz_im + ak1.im * cz_re;
		fdn = 0.0;
		if(dnu2 > 4.7170688552396617E-153)
			fdn = dnu2 * dnu2;

		ez_re = 8.0 * z.re;
		ez_im = 8.0 * z.im;
		aez = 8.0 * rt_hypotd_snf(z.re, z.im);
		s = tol / aez;
		dnu2 = rl + rl;
		if(dnu2 < 0.0)
			b_dnu2 = ceil(dnu2);
		else
			b_dnu2 = floor(dnu2);

		jl = (int)b_dnu2 + 2;
		if(z.im != 0.0)
		{
			inu = (int)fnu;
			acz = (fnu - (double)inu) * 3.1415926535897931;
			dnu2 = sin(acz);
			bk = cos(acz);
			if(z.im < 0.0)
				bk = -bk;

			if((inu & 1) != 0)
			{
				p1_re = -(-dnu2);
				bk = -bk;
			}
			else
				p1_re = -dnu2;
		}
		else
		{
			p1_re = 0.0;
			bk = 0.0;
		}

		inu = 1;
		exitg1 = false;
		while((!exitg1) && (inu <= n))
		{
			sqk = fdn - 1.0;
			b_atol = s * fabs(fdn - 1.0);
			sgn = 1.0;
			cs1_re = 1.0;
			cs1_im = 0.0;
			cs2_re = 1.0;
			cs2_im = 0.0;
			cz_re = 1.0;
			cz_im = 0.0;
			ak = 0.0;
			aa = 1.0;
			bb = aez;
			dk_re = ez_re;
			dk_im = ez_im;
			errflag = true;
			inu = 1;
			exitg2 = false;
			while((!exitg2) && (inu <= jl))
			{
				cz_re *= sqk;
				cz_im *= sqk;
				b_cz_re = cz_re;
				if(dk_im == 0.0)
				{
					if(cz_im == 0.0)
					{
						cz_re /= dk_re;
						cz_im = 0.0;
					}
					else if (cz_re == 0.0)
					{
						cz_re = 0.0;
						cz_im /= dk_re;
					}
					else
					{
						cz_re /= dk_re;
						cz_im /= dk_re;
					}
				}
				else if(dk_re == 0.0)
				{
					if(cz_re == 0.0)
					{
						cz_re = cz_im / dk_im;
						cz_im = 0.0;
					}
					else if (cz_im == 0.0)
					{
						cz_re = 0.0;
						cz_im = -(b_cz_re / dk_im);
					}
					else
					{
						cz_re = cz_im / dk_im;
						cz_im = -(b_cz_re / dk_im);
					}
				}
				else
				{
					brm = fabs(dk_re);
					acz = fabs(dk_im);
					if(brm > acz)
					{
						dnu2 = dk_im / dk_re;
						acz = dk_re + dnu2 * dk_im;
						cz_re = (cz_re + dnu2 * cz_im) / acz;
						cz_im = (cz_im - dnu2 * b_cz_re) / acz;
					}
					else if(acz == brm)
					{
						if(dk_re > 0.0)
							acz = 0.5;
						else
							acz = -0.5;

						if(dk_im > 0.0)
							dnu2 = 0.5;
						else
							dnu2 = -0.5;

						cz_re = (cz_re * acz + cz_im * dnu2) / brm;
						cz_im = (cz_im * acz - b_cz_re * dnu2) / brm;
					}
					else
					{
						dnu2 = dk_re / dk_im;
						acz = dk_im + dnu2 * dk_re;
						cz_re = (dnu2 * cz_re + cz_im) / acz;
						cz_im = (dnu2 * cz_im - b_cz_re) / acz;
					}
				}

				cs2_re += cz_re;
				cs2_im += cz_im;
				sgn = -sgn;
				cs1_re += cz_re * sgn;
				cs1_im += cz_im * sgn;
				dk_re += ez_re;
				dk_im += ez_im;
				aa = aa * fabs(sqk) / bb;
				bb += aez;
				ak += 8.0;
				sqk -= ak;
				if(aa <= b_atol)
				{
					errflag = false;
					exitg2 = true;
				}
				else
					inu++;
			}

			if(errflag)
			{
				nz = -2;
				exitg1 = true;
			}
			else
			{
				if(x + x < elim)
				{
					cz_re = -2.0 * z.re;
					cz_im = -2.0 * z.im;
					if(rtIsInf(cz_im) && rtIsInf(cz_re) && (cz_re < 0.0))
					{
						cz_re = 0.0;
						cz_im = 0.0;
					}
					else
					{
						acz = exp(cz_re / 2.0);
						cz_re = acz * (acz * cos(cz_im));
						cz_im = acz * (acz * sin(cz_im));
					}

					b_cz_re = cz_re * cs2_re - cz_im * cs2_im;
					cz_im = cz_re * cs2_im + cz_im * cs2_re;
					cs1_re += b_cz_re * p1_re - cz_im * bk;
					cs1_im += b_cz_re * bk + cz_im * p1_re;
				}

				fdn = (fdn + 8.0 * fnu) + 4.0;
				p1_re = -p1_re;
				bk = -bk;
				y->re = cs1_re * ak1.re - cs1_im * ak1.im;
				y->im = cs1_re * ak1.im + cs1_im * ak1.re;
				inu = 2;
			}
		}
	}

	return nz;
}

void cbesi(complex z, double fnu, int kode, complex *cy, int *nz, int *ierr)
{
	double AZ;
	complex zn;
	double csgn_re;
	double csgn_im;
	int inu;
	int nn;
	double dfnu;
	bool guard1 = false;
	bool guard2 = false;
	bool guard3 = false;
	int nw;
	complex cw[2];
	double b_AZ;
	int b_dfnu;
	*ierr = 0;
	AZ = rt_hypotd_snf(z.re, z.im);
	if((AZ > 1.0737418235E+9) || (fnu > 1.0737418235E+9))
	{
		printf("ierr = 4\n");
		*ierr = 4;
	}
	else
	{
		if((AZ > 32767.999992370605) || (fnu > 32767.999992370605))
		{
			printf("ierr = 3\n");
			*ierr = 3;
		}
	}

	zn = z;
	csgn_re = 1.0;
	csgn_im = 0.0;
	if(z.re < 0.0)
	{
		zn.re = -z.re;
		zn.im = -z.im;
		inu = (int)fnu;
		AZ = (fnu - (double)inu) * 3.1415926535897931;
		if(z.im < 0.0)
		{
			AZ = -AZ;
		}

		csgn_re = cos(AZ);
		csgn_im = sin(AZ);
		if((inu & 1) == 1)
		{
			csgn_re = -csgn_re;
			csgn_im = -csgn_im;
		}
	}

	cy->re = 0.0;
	cy->im = 0.0;
	*nz = 0;
	AZ = rt_hypotd_snf(zn.re, zn.im);
	nn = 1;
	dfnu = fnu;
	guard1 = false;
	guard2 = false;
	guard3 = false;
	if((AZ <= 2.0) || (AZ * AZ * 0.25 <= fnu + 1.0))
	{
		printf("cseri\n");
		nw = cseri(zn, fnu, kode, 1, cy, 2.2204460492503131E-16, 700.92179369444591, 664.87164553371019);
		//if(nw <= MIN_int32_T) {
		if(nw <= (int)(-2147483648))
		{
			//inu = MAX_int32_T;
			inu = (int)(2147483647);
		}
		else
			inu = -nw;

		if(nw < 0)
			*nz = inu;
		else
			*nz = nw;

		nn = 1 - *nz;
		if((nn == 0) || (nw >= 0))
			;
		else
		{
			dfnu = (fnu + (double)nn) - 1.0;
			guard3 = true;
		}
	}
	else
		guard3 = true;

	if(guard3)
	{
		if(AZ < 21.784271729432426)
		{
			if(dfnu <= 1.0)
			{
				printf("cmlri\n");
				nw = cmlri(zn, fnu, kode, nn, cy, 2.2204460492503131E-16);
				if(nw < 0)
				{
					if(nw == -2)
						*nz = -2;
					else
						*nz = -1;
				}
				else
					*nz = 0;
			}
			else
				guard2 = true;
		}
		else if((dfnu <= 1.0) || (!(2.0 * AZ < dfnu * dfnu)))
		{
			printf("casyi\n");
			nw = casyi(zn, fnu, kode, nn, cy, 21.784271729432426, 2.2204460492503131E-16, 700.92179369444591);
			if(nw < 0)
			{
				if(nw == -2)
					*nz = -2;
				else
					*nz = -1;
			}
			else
				*nz = 0;
		}
		else
			guard2 = true;
	}

	if(guard2)
	{
		printf("cuoik\n");
		nw = cuoik(zn, fnu, kode, 1, nn, cy, 2.2204460492503131E-16, 700.92179369444591, 664.87164553371019);
		if(nw < 0)
		{
			if(nw == -2)
				*nz = -2;
			else
				*nz = -1;
		}
		else
		{
			*nz += nw;
			nn -= nw;
			if(nn == 0)
				;
			else
			{
				dfnu = (fnu + (double)nn) - 1.0;
				if((dfnu > 85.921358647162123) || (AZ > 85.921358647162123))
				{
					if(dfnu > 86.921358647162123)
						b_dfnu = 0;
					else
						b_dfnu = (int)((85.921358647162123 - dfnu) + 1.0);

					printf("cbuni\n");
					cbuni(zn, fnu, kode, nn, cy, b_dfnu, 85.921358647162123, 2.2204460492503131E-16, 700.92179369444591, 664.87164553371019, &inu, &nw);
					if(nw < 0)
					{
						if(nw == -2)
							*nz = -2;
						else
							*nz = -1;
					}
					else
					{
						*nz += nw;
						if(inu == 0)
							;
						else
						{
							nn = inu;
							guard1 = true;
						}
					}
				}
				else
					guard1 = true;
			}
		}
	}

	if(guard1)
	{
		if(AZ > 21.784271729432426)
		{
			for(inu = 0; inu < 2; inu++)
			{
				cw[inu].re = 0.0;
				cw[inu].im = 0.0;
			}

			printf("b_cuoik\n");
			nw = b_cuoik(zn, fnu, kode, 2, 2, cw, 2.2204460492503131E-16, 700.92179369444591, 664.87164553371019);
			if(nw < 0)
			{
				*nz = nn;
				inu = 1;
				while(inu <= nn)
				{
					cy->re = 0.0;
					cy->im = 0.0;
					inu = 2;
				}
			}
			else if(nw > 0)
				*nz = -1;
			else
			{
				printf("cwrsk\n");
				nw = cwrsk(zn, fnu, kode, nn, cy, cw, 2.2204460492503131E-16, 700.92179369444591, 664.87164553371019);
				if(nw < 0)
				{
					if(nw == -2)
						*nz = -2;
					else
						*nz = -1;
				}
				else
					*nz = 0;
			}
		}
		else
		{
			printf("cmlri2\n");
			nw = cmlri(zn, fnu, kode, nn, cy, 2.2204460492503131E-16);
			if(nw < 0)
			{
				if(nw == -2)
					*nz = -2;
				else
					*nz = -1;
			}
			else
				*nz = 0;
		}
	}

	if(*nz < 0)
	{
		if(*nz == -2)
		{
			*nz = 0;
			*ierr = 5;
		}
		else
		{
			*nz = 0;
			*ierr = 2;
		}
	}
	else if((z.re >= 0.0) || (!(*nz != 1)))
		;
	else
	{
		AZ = fabs(cy->re);
		dfnu = fabs(cy->im);
		if((AZ >= dfnu) || rtIsNaN(dfnu))
		{
			printf("IsNaN detection\n");
			b_AZ = AZ;
		}
		else
			b_AZ = dfnu;

		if(b_AZ <= 1.0020841800044864E-289)
		{
			cy->re *= 4.503599627370496E+15;
			cy->im *= 4.503599627370496E+15;
			AZ = 2.2204460492503131E-16;
		}
		else
			AZ = 1.0;

		dfnu = cy->re;
		cy->re = cy->re * csgn_re - cy->im * csgn_im;
		cy->im = dfnu * csgn_im + cy->im * csgn_re;
		cy->re *= AZ;
		cy->im *= AZ;
	}
}

#if 0
void b_cuni2(const complex z, double fnu, int kode, int nin, complex y[2],
		double fnul, double tol, double elim, double alim, int
		*nlast, int *nz);
void cuni1(const complex z, double fnu, int kode, int nin, complex *y,
		double tol, double elim, double alim, int *nlast, int *nz);
void cuni2(const complex z, double fnu, int kode, int nin, complex *y,
		double tol, double elim, double alim, int *nlast, int *nz);
#endif

/* Function Definitions */
void b_cuni2(complex z, double fnu, int kode, int nin, complex y[2], double fnul, double tol, double elim, double alim, int *nlast, int *nz)
{
	int n;
	double cssr[3];
	double csrr[3];
	double bry1;
	double yy;
	complex zn;
	double zb_re;
	double zb_im;
	signed char cid_im;
	double ffnu;
	int inu;
	double c2_re;
	double c2_im;
	double zar_re;
	double zar_im;
	int in;
	static const char_complex icv1[4] = { { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } };

	double b_fnu;
	complex dai;
	complex s1;
	complex zeta1;
	complex zeta2;
	complex ai;
	complex unusedU7;
	double br;
	double bi;
	double rs1;
	double brm;
	int exitg1;
	bool goto_mw120;
	int i;
	double b_br;
	int nn;
	double b_bi;
	bool exitg2;
	double fn;
	complex phi;
	complex asum;
	complex bsum;
	bool guard1 = false;

	if(nin <= 2)
		n = nin;
	else
		n = 2;

	*nz = 0;
	*nlast = 0;
	cssr[0] = 1.0 / tol;
	cssr[1] = 1.0;
	cssr[2] = tol;
	csrr[0] = tol;
	csrr[1] = 1.0;
	csrr[2] = 1.0 / tol;
	bry1 = 2.2250738585072014E-305 / tol;
	yy = z.im;
	zn.re = z.im;
	zn.im = -z.re;
	zb_re = z.re;
	zb_im = z.im;
	cid_im = -1;

	if(fnu < 0.0)
		ffnu = ceil(fnu);
	else
		ffnu = floor(fnu);

	inu = (int)ffnu - 1;
	ffnu = 1.5707963267948966 * (fnu - ffnu);
	c2_re = cos(ffnu);
	c2_im = sin(ffnu);
	zar_re = c2_re;
	zar_im = c2_im;
	in = inu + n;
	in -= (in >> 2) << 2;
	in++;
	ffnu = c2_re;
	c2_re = c2_re * (double)icv1[in - 1].re - c2_im * (double)icv1[in - 1].im;
	c2_im = ffnu * (double)icv1[in - 1].im + c2_im * (double)icv1[in - 1].re;

	if(z.im <= 0.0)
	{
		zn.re = -z.im;
		zn.im = -z.re;
		zb_re = z.re;
		zb_im = -z.im;
		cid_im = 1;
		c2_im = -c2_im;
	}

	if (fnu >= 1.0)
		b_fnu = fnu;
	else
		b_fnu = 1.0;

	b_cunhj(zn, b_fnu, 1, tol, &dai, &s1, &zeta1, &zeta2, &ai, &unusedU7);

	if(kode == 1)
		s1.re = zeta2.re - zeta1.re;
	else
	{
		br = zb_re + zeta2.re;
		bi = zb_im + zeta2.im;

		if (bi == 0.0)
			ffnu = fnu / br;
		else if(br == 0.0)
		{
			if (fnu == 0.0)
				ffnu = 0.0 / bi;
			else
				ffnu = 0.0;
		}
		else
		{
			brm = fabs(br);
			ffnu = fabs(bi);

			if(brm > ffnu)
			{
				brm = bi / br;
				ffnu = (fnu + brm * 0.0) / (br + brm * bi);
			}
			else if(ffnu == brm)
			{
				if(br > 0.0)
					b_br = 0.5;
				else
					b_br = -0.5;

				if(bi > 0.0)
					b_bi = 0.5;
				else
					b_bi = -0.5;

				ffnu = (fnu * b_br + 0.0 * b_bi) / brm;
			}
			else
			{
				brm = br / bi;
				ffnu = brm * fnu / (bi + brm * br);
			}
		}

		s1.re = fnu * ffnu - zeta1.re;
	}

	rs1 = s1.re;
	if(fabs(s1.re) > elim)
	{
		if(s1.re > 0.0)
			*nz = -1;
		else
		{
			*nz = n;
			for(i = 1; i <= n; i++)
			{
				y[i - 1].re = 0.0;
				y[i - 1].im = 0.0;
			}
		}
	}
	else
	{
		do
		{
			exitg1 = 0;
			goto_mw120 = false;
			in = -1;
			if(2 <= n)
				nn = 2;
			else
				nn = n;

			i = 1;
			exitg2 = false;
			while((!exitg2) && (i <= nn))
			{
				fn = fnu + (double)(n - i);
				b_cunhj(zn, fn, 0, tol, &phi, &unusedU7, &zeta1, &zeta2, &asum, &bsum);
				if(kode == 1)
				{
					s1.re = zeta2.re - zeta1.re;
					s1.im = zeta2.im - zeta1.im;
				}
				else
				{
					br = zb_re + zeta2.re;
					bi = zb_im + zeta2.im;
					if(bi == 0.0)
					{
						rs1 = fn / br;
						ffnu = 0.0;
					}
					else if(br == 0.0)
					{
						if(fn == 0.0)
						{
							rs1 = 0.0 / bi;
							ffnu = 0.0;
						}
						else
						{
							rs1 = 0.0;
							ffnu = -(fn / bi);
						}
					}
					else
					{
						brm = fabs(br);
						ffnu = fabs(bi);
						if(brm > ffnu)
						{
							brm = bi / br;
							ffnu = br + brm * bi;
							rs1 = (fn + brm * 0.0) / ffnu;
							ffnu = (0.0 - brm * fn) / ffnu;
						}
						else if(ffnu == brm)
						{
							if(br > 0.0)
								br = 0.5;
							else
								br = -0.5;

							if(bi > 0.0)
								ffnu = 0.5;
							else
								ffnu = -0.5;

							rs1 = (fn * br + 0.0 * ffnu) / brm;
							ffnu = (0.0 * br - fn * ffnu) / brm;
						}
						else
						{
							brm = br / bi;
							ffnu = bi + brm * br;
							rs1 = brm * fn / ffnu;
							ffnu = (brm * 0.0 - fn) / ffnu;
						}
					}

					s1.re = (fn * rs1 - zeta1.re) + fabs(yy) * 0.0;
					s1.im = (fn * ffnu - zeta1.im) + fabs(yy);
				}

				rs1 = s1.re;
				if(fabs(s1.re) > elim)
				{
					goto_mw120 = true;
					exitg2 = true;
				}
				else
				{
					if(i == 1)
						in = 1;

					guard1 = false;
					if(fabs(s1.re) >= alim)
					{
						rs1 = ((s1.re + log(rt_hypotd_snf(phi.re, phi.im))) - 0.25 * log(rt_hypotd_snf(unusedU7.re, unusedU7.im))) - 1.2655121234846454;
						if(fabs(rs1) > elim)
						{
							goto_mw120 = true;
							exitg2 = true;
						}
						else
						{
							if(i == 1)
								in = 0;

							if((rs1 >= 0.0) && (i == 1))
								in = 2;

							guard1 = true;
						}
					}
					else
						guard1 = true;

					if(guard1)
					{
						ai = cairy(unusedU7, 0, 2);
						dai = cairy(unusedU7, 1, 2);
						ffnu = (ai.re * asum.re - ai.im * asum.im) + (dai.re * bsum.re - dai.im * bsum.im);
						brm = (ai.re * asum.im + ai.im * asum.re) + (dai.re * bsum.im + dai.im * bsum.re);
						bi = ffnu * phi.re - brm * phi.im;
						brm = ffnu * phi.im + brm * phi.re;
						ffnu = exp(s1.re) * cssr[in] * cos(s1.im);
						br = exp(s1.re) * cssr[in] * sin(s1.im);
						dai.re = bi * ffnu - brm * br;
						dai.im = bi * br + brm * ffnu;
						if((in + 1 == 1) && (cuchk(dai, bry1, tol) != 0))
						{
							goto_mw120 = true;
							exitg2 = true;
						}
						else
						{
							if(yy <= 0.0)
								dai.im = -dai.im;

							y[n - i].re = csrr[in] * (dai.re * c2_re - dai.im * c2_im);
							y[n - i].im = csrr[in] * (dai.re * c2_im + dai.im * c2_re);
							ffnu = c2_re;
							c2_re = c2_re * 0.0 - c2_im * (double)cid_im;
							c2_im = ffnu * (double)cid_im + c2_im * 0.0;
							i++;
						}
					}
				}
			}

			if(!goto_mw120)
				exitg1 = 1;
			else if(rs1 > 0.0)
			{
				*nz = -1;
				exitg1 = 1;
			}
			else
			{
				y[n - 1].re = 0.0;
				y[n - 1].im = 0.0;
				(*nz)++;
				if(n - 1 == 0)
					exitg1 = 1;
				else
				{
					in = b_cuoik(z, fnu, kode, 1, 1, y, tol, elim, alim);
					if(in < 0)
					{
						*nz = -1;
						exitg1 = 1;
					}
					else
					{
						n = 1 - in;
						*nz += in;
						if(n == 0)
							exitg1 = 1;
						else if((fnu + (double)n) - 1.0 < fnul)
						{
							*nlast = n;
							exitg1 = 1;
						}
						else
						{
							in = inu + n;
							nn = icv1[in - ((in >> 2) << 2)].re;
							in = icv1[in - ((in >> 2) << 2)].im;
							c2_re = zar_re * (double)nn - zar_im * (double)in;
							c2_im = zar_re * (double)in + zar_im * (double)nn;
							if(yy <= 0.0)
								c2_im = -c2_im;
						}
					}
				}
			}
		}
		while(exitg1 == 0);
	}
}

void cuni1(complex z, double fnu, int kode, int nin, complex *y, double tol, double elim, double alim, int *nlast, int *nz)
{
	int n;
	double cssr[3];
	double csrr[3];
	double bry1;
	double fn;
	int iflag;
	complex cwrk[16];
	complex s1;
	complex zeta1;
	complex zeta2;
	complex unusedU1;
	double br;
	double bi;
	double rs1;
	double brm;
	bool goto_mw110;
	bool exitg1;
	double b_br;
	complex summ;
	double b_bi;
	bool guard1 = false;
	double sgnbr;

	if(nin <= 1)
		n = nin;
	else
		n = 1;

	*nz = 0;
	*nlast = 0;
	cssr[0] = 1.0 / tol;
	cssr[1] = 1.0;
	cssr[2] = tol;
	csrr[0] = tol;
	csrr[1] = 1.0;
	csrr[2] = 1.0 / tol;
	bry1 = 2.2250738585072014E-305 / tol;

	if(fnu >= 1.0)
		fn = fnu;
	else
		fn = 1.0;

	iflag = 0;
	memset(&cwrk[0], 0, sizeof(complex) << 4);
	b_cunik(z, fn, 1, 1, tol, &iflag, cwrk, &s1, &zeta1, &zeta2, &unusedU1);

	if(kode == 1)
		s1.re = zeta2.re - zeta1.re;
	else
	{
		br = z.re + zeta2.re;
		bi = z.im + zeta2.im;

		if(bi == 0.0)
			br = fn / br;
		else if(br == 0.0)
		{
			if(fn == 0.0)
				br = 0.0 / bi;
			else
				br = 0.0;
		}
		else
		{
			brm = fabs(br);
			rs1 = fabs(bi);

			if(brm > rs1)
			{
				brm = bi / br;
				br = (fn + brm * 0.0) / (br + brm * bi);
			}
			else if(rs1 == brm)
			{
				if(br > 0.0)
					b_br = 0.5;
				else
					b_br = -0.5;

				if(bi > 0.0)
					b_bi = 0.5;
				else
					b_bi = -0.5;

				br = (fn * b_br + 0.0 * b_bi) / brm;
			}
			else
			{
				brm = br / bi;
				br = brm * fn / (bi + brm * br);
			}
		}

		s1.re = fn * br - zeta1.re;
	}

	rs1 = s1.re;
	if(fabs(s1.re) > elim)
	{
		if(s1.re > 0.0)
			*nz = -1;
		else
		{
			*nz = n;
			iflag = 1;
			while(iflag <= n)
			{
				y->re = 0.0;
				y->im = 0.0;
				iflag = 2;
			}
		}
	}
	else
	{
		goto_mw110 = false;
		iflag = 1;
		exitg1 = false;
		while((!exitg1) && (iflag <= n))
		{
			fn = fnu;
			iflag = 0;
			b_cunik(z, fn, 1, 0, tol, &iflag, cwrk, &unusedU1, &zeta1, &zeta2, &summ);
			if(kode == 1)
			{
				s1.re = zeta2.re - zeta1.re;
				s1.im = zeta2.im - zeta1.im;
			}
			else
			{
				br = z.re + zeta2.re;
				bi = z.im + zeta2.im;
				if(bi == 0.0)
				{
					br = fn / br;
					rs1 = 0.0;
				}
				else if(br == 0.0)
				{
					if(fn == 0.0)
					{
						br = 0.0 / bi;
						rs1 = 0.0;
					}
					else
					{
						br = 0.0;
						rs1 = -(fn / bi);
					}
				}
				else
				{
					brm = fabs(br);
					rs1 = fabs(bi);
					if(brm > rs1)
					{
						brm = bi / br;
						rs1 = br + brm * bi;
						br = (fn + brm * 0.0) / rs1;
						rs1 = (0.0 - brm * fn) / rs1;
					}
					else if(rs1 == brm)
					{
						if (br > 0.0)
							sgnbr = 0.5;
						else
							sgnbr = -0.5;

						if (bi > 0.0)
							rs1 = 0.5;
						else
							rs1 = -0.5;

						br = (fn * sgnbr + 0.0 * rs1) / brm;
						rs1 = (0.0 * sgnbr - fn * rs1) / brm;
					}
					else
					{
						brm = br / bi;
						rs1 = bi + brm * br;
						br = brm * fn / rs1;
						rs1 = (brm * 0.0 - fn) / rs1;
					}
				}

				s1.re = fn * br - zeta1.re;
				s1.im = (fn * rs1 - zeta1.im) + z.im;
			}

			rs1 = s1.re;
			if(fabs(s1.re) > elim)
			{
				goto_mw110 = true;
				exitg1 = true;
			}
			else
			{
				iflag = 1;
				guard1 = false;
				if(fabs(s1.re) >= alim)
				{
					rs1 = s1.re + log(rt_hypotd_snf(unusedU1.re, unusedU1.im));
					if(fabs(rs1) > elim)
					{
						goto_mw110 = true;
						exitg1 = true;
					}
					else
					{
						iflag = 0;
						if(rs1 >= 0.0)
							iflag = 2;

						guard1 = true;
					}
				}
				else
					guard1 = true;

				if(guard1)
				{
					br = unusedU1.re * summ.re - unusedU1.im * summ.im;
					brm = unusedU1.re * summ.im + unusedU1.im * summ.re;
					sgnbr = exp(s1.re) * cssr[iflag] * cos(s1.im);
					bi = exp(s1.re) * cssr[iflag] * sin(s1.im);
					s1.re = br * sgnbr - brm * bi;
					s1.im = br * bi + brm * sgnbr;
					if((iflag + 1 == 1) && (cuchk(s1, bry1, tol) != 0))
					{
						goto_mw110 = true;
						exitg1 = true;
					}
					else
					{
						y->re = csrr[iflag] * s1.re;
						y->im = csrr[iflag] * s1.im;
						iflag = 2;
					}
				}
			}
		}

		if(!goto_mw110)
			;
		else if (rs1 > 0.0)
			*nz = -1;
		else
		{
			y->re = 0.0;
			y->im = 0.0;
			*nz = 1;
		}
	}
}

void cuni2(complex z, double fnu, int kode, int nin, complex *y, double tol, double elim, double alim, int *nlast, int *nz)
{
	int n;
	double cssr[3];
	double csrr[3];
	double bry1;
	double yy;
	complex zn;
	double zb_re;
	double zb_im;
	signed char cid_im;
	double ffnu;
	double ang;
	int in;
	double re;
	double c2_re;
	static const char_complex icv0[4] = { { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } };

	double c2_im;
	double b_fnu;
	complex dai;
	complex s1;
	complex zeta1;
	complex zeta2;
	complex ai;
	complex unusedU7;
	double bi;
	double rs1;
	bool goto_mw120;
	bool exitg1;
	double b_ang;
	double fn;
	complex phi;
	complex asum;
	complex bsum;
	double b_bi;
	bool guard1 = false;

	if(nin <= 1)
		n = nin;
	else
		n = 1;

	*nz = 0;
	*nlast = 0;
	cssr[0] = 1.0 / tol;
	cssr[1] = 1.0;
	cssr[2] = tol;
	csrr[0] = tol;
	csrr[1] = 1.0;
	csrr[2] = 1.0 / tol;
	bry1 = 2.2250738585072014E-305 / tol;
	yy = z.im;
	zn.re = z.im;
	zn.im = -z.re;
	zb_re = z.re;
	zb_im = z.im;
	cid_im = -1;

	if(fnu < 0.0)
		ffnu = ceil(fnu);
	else
		ffnu = floor(fnu);

	ang = 1.5707963267948966 * (fnu - ffnu);
	in = ((int)ffnu + n) - 1;
	in -= (in >> 2) << 2;
	re = cos(ang);
	ffnu = sin(ang);
	c2_re = re * (double)icv0[in].re - ffnu * (double)icv0[in].im;
	c2_im = re * (double)icv0[in].im + ffnu * (double)icv0[in].re;

	if(z.im <= 0.0)
	{
		zn.re = -z.im;
		zn.im = -z.re;
		zb_re = z.re;
		zb_im = -z.im;
		cid_im = 1;
		c2_im = -c2_im;
	}

	if(fnu >= 1.0)
		b_fnu = fnu;
	else
		b_fnu = 1.0;

	b_cunhj(zn, b_fnu, 1, tol, &dai, &s1, &zeta1, &zeta2, &ai, &unusedU7);

	if(kode == 1)
		s1.re = zeta2.re - zeta1.re;
	else
	{
		ang = zb_re + zeta2.re;
		bi = zb_im + zeta2.im;

		if(bi == 0.0)
			ffnu = fnu / ang;
		else if(ang == 0.0)
		{
			if(fnu == 0.0)
				ffnu = 0.0 / bi;
			else
				ffnu = 0.0;
		}
		else
		{
			re = fabs(ang);
			ffnu = fabs(bi);

			if(re > ffnu)
			{
				re = bi / ang;
				ffnu = (fnu + re * 0.0) / (ang + re * bi);
			}
			else if(ffnu == re)
			{
				if(ang > 0.0)
					b_ang = 0.5;
				else
					b_ang = -0.5;

				if(bi > 0.0)
					b_bi = 0.5;
				else
					b_bi = -0.5;

				ffnu = (fnu * b_ang + 0.0 * b_bi) / re;
			}
			else
			{
				re = ang / bi;
				ffnu = re * fnu / (bi + re * ang);
			}
		}

		s1.re = fnu * ffnu - zeta1.re;
	}

	rs1 = s1.re;

	if(fabs(s1.re) > elim)
	{
		if (s1.re > 0.0)
			*nz = -1;
		else
		{
			*nz = n;
			in = 1;

			while(in <= n)
			{
				y->re = 0.0;
				y->im = 0.0;
				in = 2;
			}
		}
	}
	else
	{
		goto_mw120 = false;
		in = 1;
		exitg1 = false;

		while((!exitg1) && (in <= n))
		{
			fn = fnu;
			b_cunhj(zn, fn, 0, tol, &phi, &unusedU7, &zeta1, &zeta2, &asum, &bsum);

			if(kode == 1)
			{
				s1.re = zeta2.re - zeta1.re;
				s1.im = zeta2.im - zeta1.im;
			}
			else
			{
				ang = zb_re + zeta2.re;
				bi = zb_im + zeta2.im;

				if(bi == 0.0)
				{
					rs1 = fn / ang;
					ffnu = 0.0;
				}
				else if(ang == 0.0)
				{
					if(fn == 0.0)
					{
						rs1 = 0.0 / bi;
						ffnu = 0.0;
					}
					else
					{
						rs1 = 0.0;
						ffnu = -(fn / bi);
					}
				}
				else
				{
					re = fabs(ang);
					ffnu = fabs(bi);

					if(re > ffnu)
					{
						re = bi / ang;
						ffnu = ang + re * bi;
						rs1 = (fn + re * 0.0) / ffnu;
						ffnu = (0.0 - re * fn) / ffnu;
					}
					else if(ffnu == re)
					{
						if (ang > 0.0)
							ang = 0.5;
						else
							ang = -0.5;

						if (bi > 0.0)
							ffnu = 0.5;
						else
							ffnu = -0.5;

						rs1 = (fn * ang + 0.0 * ffnu) / re;
						ffnu = (0.0 * ang - fn * ffnu) / re;
					}
					else
					{
						re = ang / bi;
						ffnu = bi + re * ang;
						rs1 = re * fn / ffnu;
						ffnu = (re * 0.0 - fn) / ffnu;
					}
				}

				s1.re = (fn * rs1 - zeta1.re) + fabs(yy) * 0.0;
				s1.im = (fn * ffnu - zeta1.im) + fabs(yy);
			}

			rs1 = s1.re;
			if(fabs(s1.re) > elim)
			{
				goto_mw120 = true;
				exitg1 = true;
			}
			else
			{
				in = 1;
				guard1 = false;

				if(fabs(s1.re) >= alim)
				{
					rs1 = ((s1.re + log(rt_hypotd_snf(phi.re, phi.im))) - 0.25 * log(rt_hypotd_snf(unusedU7.re, unusedU7.im))) - 1.2655121234846454;

					if(fabs(rs1) > elim)
					{
						goto_mw120 = true;
						exitg1 = true;
					}
					else
					{
						in = 0;

						if(rs1 >= 0.0)
							in = 2;

						guard1 = true;
					}
				}
				else
					guard1 = true;

				if(guard1)
				{
					ai = cairy(unusedU7, 0, 2);
					dai = cairy(unusedU7, 1, 2);
					ffnu = (ai.re * asum.re - ai.im * asum.im) + (dai.re * bsum.re - dai.im * bsum.im);
					ang = (ai.re * asum.im + ai.im * asum.re) + (dai.re * bsum.im + dai.im * bsum.re);
					bi = ffnu * phi.re - ang * phi.im;
					ang = ffnu * phi.im + ang * phi.re;
					re = exp(s1.re) * cssr[in] * cos(s1.im);
					ffnu = exp(s1.re) * cssr[in] * sin(s1.im);
					dai.re = bi * re - ang * ffnu;
					dai.im = bi * ffnu + ang * re;

					if((in + 1 == 1) && (cuchk(dai, bry1, tol) != 0))
					{
						goto_mw120 = true;
						exitg1 = true;
					}
					else
					{
						if(yy <= 0.0)
							dai.im = -dai.im;

						y->re = csrr[in] * (dai.re * c2_re - dai.im * c2_im);
						y->im = csrr[in] * (dai.re * c2_im + dai.im * c2_re);
						ffnu = c2_re;
						c2_re = c2_re * 0.0 - c2_im * (double)cid_im;
						c2_im = ffnu * (double)cid_im + c2_im * 0.0;
						in = 2;
					}
				}
			}
		}

		if(!goto_mw120)
			;
		else if(rs1 > 0.0)
			*nz = -1;
		else
		{
			y->re = 0.0;
			y->im = 0.0;
			*nz = 1;
		}
	}
}

void cbuni(complex z, double fnu, int kode, int nin, complex *y, int nui, double fnul, double tol, double elim, double alim, int *nlast, int *nz)
{
	int n;
	int iform;
	double fnui;
	double dfnu;
	int nw;
	double gnu;
	complex cy[2];
	int nd;
	double cssr[3];
	double rs1;
	double bry0;
	double csrr[3];
	double bry1;
	int iflag;
	double fn;
	double dscl;
	complex cwrk[16];
	complex s2;
	complex s1;
	complex zeta2;
	complex rz;
	double dscr;
	double brm;
	int exitg1;
	int nn;
	double b_bry0;
	bool goto_mw110;
	double b_rs1;
	bool exitg2;
	int unusedU3;
	bool guard1 = false;
	double c_rs1;

	if(nin <= 1)
		n = nin;
	else
		n = 1;

	*nz = 0;
	iform = 1;

	if(fabs(z.im) > fabs(z.re) * 1.7321)
		iform = 2;

	if(nui == 0)
	{
		if(iform == 2)
			cuni2(z, fnu, kode, n, y, tol, elim, alim, nlast, &nw);
		else
			cuni1(z, fnu, kode, n, y, tol, elim, alim, nlast, &nw);

		if(nw < 0)
		{
			*nz = -1;

			if(nw == -2)
				*nz = -2;
		}
		else
			*nz = nw;
	}
	else
	{
		fnui = nui;
		dfnu = (fnu + (double)n) - 1.0;
		gnu = dfnu + (double)nui;

		if(iform == 2)
		{
			for(iform = 0; iform < 2; iform++)
			{
				cy[iform].re = 0.0;
				cy[iform].im = 0.0;
			}

			b_cuni2(z, gnu, kode, 2, cy, fnul, tol, elim, alim, nlast, &nw);
		}
		else
		{
			for(iform = 0; iform < 2; iform++)
			{
				cy[iform].re = 0.0;
				cy[iform].im = 0.0;
			}

			nw = 0;
			nd = 2;
			*nlast = 0;
			cssr[0] = 1.0 / tol;
			cssr[1] = 1.0;
			cssr[2] = tol;
			csrr[0] = tol;
			csrr[1] = 1.0;
			csrr[2] = 1.0 / tol;
			bry1 = 2.2250738585072014E-305 / tol;

			if(gnu >= 1.0)
				fn = gnu;
			else
				fn = 1.0;

			iform = 0;
			memset(&cwrk[0], 0, sizeof(complex) << 4);
			b_cunik(z, fn, 1, 1, tol, &iform, cwrk, &s2, &s1, &zeta2, &rz);

			if(kode == 1)
				s1.re = zeta2.re - s1.re;
			else
			{
				bry0 = z.re + zeta2.re;
				rs1 = z.im + zeta2.im;

				if(rs1 == 0.0)
				{
					dscr = fn / bry0;
				}
				else if(bry0 == 0.0)
				{
					if(fn == 0.0)
						dscr = 0.0 / rs1;
					else
						dscr = 0.0;
				}
				else
				{
					brm = fabs(bry0);
					dscl = fabs(rs1);

					if(brm > dscl)
					{
						dscl = rs1 / bry0;
						dscr = (fn + dscl * 0.0) / (bry0 + dscl * rs1);
					}
					else if(dscl == brm)
					{
						if(bry0 > 0.0)
							b_bry0 = 0.5;
						else
							b_bry0 = -0.5;

						if(rs1 > 0.0)
							b_rs1 = 0.5;
						else
							b_rs1 = -0.5;

						dscr = (fn * b_bry0 + 0.0 * b_rs1) / brm;
					}
					else
					{
						dscl = bry0 / rs1;
						dscr = dscl * fn / (rs1 + dscl * bry0);
					}
				}

				s1.re = fn * dscr - s1.re;
			}

			rs1 = s1.re;

			if(fabs(s1.re) > elim)
			{
				if(s1.re > 0.0)
					nw = -1;
				else
				{
					nw = 2;

					for(iform = 0; iform < 2; iform++)
					{
						cy[iform].re = 0.0;
						cy[iform].im = 0.0;
					}
				}
			}
			else
			{
				do
				{
					exitg1 = 0;
					iflag = -1;

					if(2 <= nd)
						nn = 2;
					else
						nn = nd;

					goto_mw110 = false;
					iform = 1;
					exitg2 = false;

					while((!exitg2) && (iform <= nn))
					{
						fn = gnu + (double)(nd - iform);
						unusedU3 = 0;
						b_cunik(z, fn, 1, 0, tol, &unusedU3, cwrk, &s2, &s1, &zeta2, &rz);

						if(kode == 1)
						{
							s1.re = zeta2.re - s1.re;
							s1.im = zeta2.im - s1.im;
						}
						else
						{
							bry0 = z.re + zeta2.re;
							rs1 = z.im + zeta2.im;

							if(rs1 == 0.0)
							{
								dscr = fn / bry0;
								rs1 = 0.0;
							}
							else if(bry0 == 0.0)
							{
								if(fn == 0.0)
								{
									dscr = 0.0 / rs1;
									rs1 = 0.0;
								}
								else
								{
									dscr = 0.0;
									rs1 = -(fn / rs1);
								}
							}
							else
							{
								brm = fabs(bry0);
								dscl = fabs(rs1);

								if(brm > dscl)
								{
									dscl = rs1 / bry0;
									bry0 += dscl * rs1;
									dscr = (fn + dscl * 0.0) / bry0;
									rs1 = (0.0 - dscl * fn) / bry0;
								}
								else if(dscl == brm)
								{
									if(bry0 > 0.0)
										dscl = 0.5;
									else
										dscl = -0.5;

									if(rs1 > 0.0)
										bry0 = 0.5;
									else
										bry0 = -0.5;

									dscr = (fn * dscl + 0.0 * bry0) / brm;
									rs1 = (0.0 * dscl - fn * bry0) / brm;
								}
								else
								{
									dscl = bry0 / rs1;
									bry0 = rs1 + dscl * bry0;
									dscr = dscl * fn / bry0;
									rs1 = (dscl * 0.0 - fn) / bry0;
								}
							}

							s1.re = fn * dscr - s1.re;
							s1.im = fn * rs1 - s1.im;
							s1.im += z.im;
						}

						rs1 = s1.re;

						if(fabs(s1.re) > elim)
						{
							goto_mw110 = true;
							exitg2 = true;
						}
						else
						{
							if(iform == 1)
								iflag = 1;

							guard1 = false;

							if(fabs(s1.re) >= alim)
							{
								rs1 = s1.re + log(rt_hypotd_snf(s2.re, s2.im));

								if(fabs(rs1) > elim)
								{
									goto_mw110 = true;
									exitg2 = true;
								}
								else
								{
									if(iform == 1)
										iflag = 0;

									if((rs1 >= 0.0) && (iform == 1))
										iflag = 2;

									guard1 = true;
								}
							}
							else
								guard1 = true;

							if(guard1)
							{
								brm = s2.re * rz.re - s2.im * rz.im;
								bry0 = s2.re * rz.im + s2.im * rz.re;
								dscl = exp(s1.re) * cssr[iflag] * cos(s1.im);
								dscr = exp(s1.re) * cssr[iflag] * sin(s1.im);
								s2.re = brm * dscl - bry0 * dscr;
								s2.im = brm * dscr + bry0 * dscl;

								if((iflag + 1 == 1) && (cuchk(s2, bry1, tol) != 0))
								{
									goto_mw110 = true;
									exitg2 = true;
								}
								else
								{
									cy[nd - iform].re = csrr[iflag] * s2.re;
									cy[nd - iform].im = csrr[iflag] * s2.im;
									iform++;
								}
							}
						}
					}

					if(!goto_mw110)
						exitg1 = 1;
					else if(rs1 > 0.0)
					{
						nw = -1;
						exitg1 = 1;
					}
					else
					{
						cy[nd - 1].re = 0.0;
						cy[nd - 1].im = 0.0;
						nw++;

						if(nd - 1 == 0)
							exitg1 = 1;
						else
						{
							iform = b_cuoik(z, gnu, kode, 1, 1, cy, tol, elim, alim);

							if(iform < 0)
							{
								nw = -1;
								exitg1 = 1;
							}
							else
							{
								nd = 1 - iform;
								nw += iform;

								if(nd == 0)
									exitg1 = 1;
								else
								{
									if(!((gnu + (double)nd) - 1.0 >= fnul))
									{
										*nlast = nd;
										exitg1 = 1;
									}
								}
							}
						}
					}
				}
				while(exitg1 == 0);
			}
		}

		if(nw < 0)
		{
			*nz = -1;

			if(nw == -2)
				*nz = -2;
		}
		else if(nw != 0)
			*nlast = n;
		else
		{
			rs1 = rt_hypotd_snf(cy[0].re, cy[0].im);
			bry0 = 2.2250738585072014E-305 / tol;
			bry1 = 1.0 / bry0;
			cssr[0] = bry0;
			cssr[1] = bry1;
			cssr[2] = bry1;
			iflag = 2;
			fn = 1.0;
			dscl = 1.0;

			if(rs1 > bry0)
			{
				if(rs1 >= bry1)
				{
					iflag = 3;
					fn = tol;
					dscl = tol;
				}
			}
			else
			{
				iflag = 1;
				bry1 = bry0;
				fn = 1.0 / tol;
				dscl = fn;
			}

			dscr = 1.0 / fn;
			s1.re = dscl * cy[1].re;
			s1.im = dscl * cy[1].im;
			s2.re = dscl * cy[0].re;
			s2.im = dscl * cy[0].im;

			if(z.im == 0.0)
			{
				rz.re = 2.0 / z.re;
				rz.im = 0.0;
			}
			else if(z.re == 0.0)
			{
				rz.re = 0.0;
				rz.im = -(2.0 / z.im);
			}
			else
			{
				brm = fabs(z.re);
				dscl = fabs(z.im);

				if(brm > dscl)
				{
					dscl = z.im / z.re;
					bry0 = z.re + dscl * z.im;
					rz.re = (2.0 + dscl * 0.0) / bry0;
					rz.im = (0.0 - dscl * 2.0) / bry0;
				}
				else if(dscl == brm)
				{
					if(z.re > 0.0)
						dscl = 0.5;
					else
						dscl = -0.5;

					if(z.im > 0.0)
						bry0 = 0.5;
					else
						bry0 = -0.5;

					rz.re = (2.0 * dscl + 0.0 * bry0) / brm;
					rz.im = (0.0 * dscl - 2.0 * bry0) / brm;
				}
				else
				{
					dscl = z.re / z.im;
					bry0 = z.im + dscl * z.re;
					rz.re = dscl * 2.0 / bry0;
					rz.im = (dscl * 0.0 - 2.0) / bry0;
				}
			}

			for(iform = 1; iform <= nui; iform++)
			{
				*y = s2;
				brm = s2.re;
				s2.re = s2.re * rz.re - s2.im * rz.im;
				s2.im = brm * rz.im + s2.im * rz.re;
				s2.re *= dfnu + fnui;
				s2.im *= dfnu + fnui;
				s2.re += s1.re;
				s2.im += s1.im;
				s1 = *y;
				fnui--;

				if(!(iflag >= 3))
				{
					y->re = dscr * s2.re;
					y->im = dscr * s2.im;
					rs1 = fabs(y->re);
					bry0 = fabs(y->im);

					if((rs1 >= bry0) || rtIsNaN(bry0))
						c_rs1 = rs1;
					else
						c_rs1 = bry0;

					if(!(c_rs1 <= bry1))
					{
						iflag++;
						bry1 = cssr[iflag - 1];
						s1.re *= dscr;
						s1.im *= dscr;
						fn *= tol;
						dscr = 1.0 / fn;
						s1.re *= fn;
						s1.im *= fn;
						s2.re = fn * y->re;
						s2.im = fn * y->im;
					}
				}
			}

			y->re = dscr * s2.re;
			y->im = dscr * s2.im;
		}
	}
}

int cwrsk(complex zr, double fnu, int kode, int nin, complex *y, complex cw[2], double tol, double elim, double alim)
{
	int nz;
	int n;
	complex b_cw[2];
	int magz;
	int itime;
	long long i3;
	double p1_re;
	double p1_im;
	double flam;
	int id;
	double test;
	double pt_re;
	int k;
	double rz_re;
	double brm;
	double rz_im;
	double p2_re;
	double p2_im;
	double b_magz;
	double t1_re;
	double c_magz;
	double t1_im;
	double ap2;
	double test1;
	double ap1;
	double pt_im;

	if(nin <= 1)
		n = nin;
	else
		n = 1;

	nz = 0;

	for(magz = 0; magz < 2; magz++)
		b_cw[magz] = cw[magz];

	itime = b_cbknu(zr, fnu, kode, 2, b_cw, tol, elim, alim);

	if(itime != 0)
	{
		nz = -1;

		if(itime == -2)
			nz = -2;
	}
	else
	{
		if(n < 1)
			;
		else
		{
			n = (int)fnu;
			magz = (int)rt_hypotd_snf(zr.re, zr.im);
			i3 = magz + 1LL;

			if(i3 > 2147483647LL)
				i3 = 2147483647LL;
			else
			{
				if (i3 < -2147483648LL)
					i3 = -2147483648LL;
			}

			id = (int)i3;

			if(n > id)
				id = 0;
			else
				id -= n;

			k = 1;

			if(zr.im == 0.0)
			{
				rz_re = 2.0 / zr.re;
				rz_im = 0.0;
			}
			else if(zr.re == 0.0)
			{
				rz_re = 0.0;
				rz_im = -(2.0 / zr.im);
			}
			else
			{
				brm = fabs(zr.re);
				flam = fabs(zr.im);

				if(brm > flam)
				{
					test = zr.im / zr.re;
					flam = zr.re + test * zr.im;
					rz_re = (2.0 + test * 0.0) / flam;
					rz_im = (0.0 - test * 2.0) / flam;
				}
				else if(flam == brm)
				{
					if(zr.re > 0.0)
						test = 0.5;
					else
						test = -0.5;

					if(zr.im > 0.0)
						flam = 0.5;
					else
						flam = -0.5;

					rz_re = (2.0 * test + 0.0 * flam) / brm;
					rz_im = (0.0 * test - 2.0 * flam) / brm;
				}
				else
				{
					test = zr.re / zr.im;
					flam = zr.im + test * zr.re;
					rz_re = test * 2.0 / flam;
					rz_im = (test * 0.0 - 2.0) / flam;
				}
			}

			if((double)magz + 1.0 >= n)
				b_magz = (double)magz + 1.0;
			else
				b_magz = n;

			t1_re = b_magz * rz_re;

			if ((double)magz + 1.0 >= n)
				c_magz = (double)magz + 1.0;
			else
				c_magz = n;

			t1_im = c_magz * rz_im;
			p2_re = -t1_re;
			p2_im = -t1_im;
			t1_re += rz_re;
			t1_im += rz_im;
			ap2 = rt_hypotd_snf(p2_re, p2_im);
			test1 = sqrt((ap2 + ap2) / tol);
			test = test1;
			p1_re = 1.0;
			p1_im = 0.0;
			itime = 1;

			while(itime <= 2)
			{
				k++;
				ap1 = ap2;
				pt_re = p2_re;
				pt_im = p2_im;
				brm = p2_re;
				p2_re = p1_re - (t1_re * p2_re - t1_im * p2_im);
				p2_im = p1_im - (t1_re * p2_im + t1_im * brm);
				p1_re = pt_re;
				p1_im = pt_im;
				t1_re += rz_re;
				t1_im += rz_im;
				ap2 = rt_hypotd_snf(p2_re, p2_im);

				if(!(ap1 <= test))
				{
					if(itime == 1)
					{
						test = rt_hypotd_snf(t1_re, t1_im) * 0.5;
						flam = test + sqrt(test * test - 1.0);
						brm = ap2 / ap1;

						if((brm <= flam) || rtIsNaN(flam))
							flam = brm;

						test = test1 * sqrt(flam / (flam * flam - 1.0));
					}

					itime++;
				}
			}

			itime = (k + id) + 1;
			t1_re = itime;
			p1_re = 1.0 / ap2;
			p1_im = 0.0;
			p2_re = 0.0;
			p2_im = 0.0;

			for(magz = 1; magz <= itime; magz++)
			{
				pt_re = p1_re;
				pt_im = p1_im;
				flam = ((fnu + 1.0) - 1.0) + t1_re;
				ap2 = rz_re * flam - rz_im * 0.0;
				test = rz_re * 0.0 + rz_im * flam;
				brm = ap2 * p1_im + test * p1_re;
				p1_re = (ap2 * p1_re - test * p1_im) + p2_re;
				p1_im = brm + p2_im;
				p2_re = pt_re;
				p2_im = pt_im;
				t1_re--;
			}

			if((p1_re == 0.0) && (p1_im == 0.0))
			{
				p1_re = tol;
				p1_im = tol;
			}

			if(p1_im == 0.0)
			{
				if(p2_im == 0.0)
				{
					y->re = p2_re / p1_re;
					y->im = 0.0;
				}
				else if(p2_re == 0.0)
				{
					y->re = 0.0;
					y->im = p2_im / p1_re;
				}
				else
				{
					y->re = p2_re / p1_re;
					y->im = p2_im / p1_re;
				}
			}
			else if(p1_re == 0.0)
			{
				if(p2_re == 0.0)
				{
					y->re = p2_im / p1_im;
					y->im = 0.0;
				}
				else if(p2_im == 0.0)
				{
					y->re = 0.0;
					y->im = -(p2_re / p1_im);
				}
				else
				{
					y->re = p2_im / p1_im;
					y->im = -(p2_re / p1_im);
				}
			}
			else
			{
				brm = fabs(p1_re);
				flam = fabs(p1_im);

				if(brm > flam)
				{
					test = p1_im / p1_re;
					flam = p1_re + test * p1_im;
					y->re = (p2_re + test * p2_im) / flam;
					y->im = (p2_im - test * p2_re) / flam;
				}
				else if(flam == brm)
				{
					if(p1_re > 0.0)
						test = 0.5;
					else
						test = -0.5;

					if(p1_im > 0.0)
						flam = 0.5;
					else
						flam = -0.5;

					y->re = (p2_re * test + p2_im * flam) / brm;
					y->im = (p2_im * test - p2_re * flam) / brm;
				}
				else
				{
					test = p1_re / p1_im;
					flam = p1_im + test * p1_re;
					y->re = (test * p2_re + p2_im) / flam;
					y->im = (test * p2_im - p2_re) / flam;
				}
			}
		}

		if(kode == 1)
		{
			p1_re = 1.0;
			p1_im = 0.0;
		}
		else
		{
			p1_re = cos(zr.im);
			p1_im = sin(zr.im);
		}

		flam = rt_hypotd_snf(b_cw[1].re, b_cw[1].im);
		test = 2.2250738585072014E-305 / tol;

		if (flam > test)
		{
			test = 1.0 / test;

			if(flam >= test)
				pt_re = tol;
			else
				pt_re = 1.0;
		}
		else
			pt_re = 1.0 / tol;

		flam = b_cw[0].re * pt_re - b_cw[0].im * 0.0;
		test = b_cw[0].re * 0.0 + b_cw[0].im * pt_re;
		brm = (y->re * flam - y->im * test) + (b_cw[1].re * pt_re - b_cw[1].im * 0.0);
		flam = (y->re * test + y->im * flam) + (b_cw[1].re * 0.0 + b_cw[1].im * pt_re);
		p2_re = brm * zr.re - flam * zr.im;
		p2_im = brm * zr.im + flam * zr.re;
		flam = 1.0 / rt_hypotd_snf(p2_re, p2_im);
		p2_im = -p2_im;
		p2_re *= flam;
		p2_im *= flam;
		p1_re *= flam;
		p1_im *= flam;
		flam = p1_re;
		p1_re = p1_re * p2_re - p1_im * p2_im;
		p1_im = flam * p2_im + p1_im * p2_re;
		y->re = p1_re * pt_re - p1_im * 0.0;
		y->im = p1_re * 0.0 + p1_im * pt_re;
	}

	return nz;
}

int b_ckscl(complex zr, int nin, complex y[2], double ascle, double tol, double elim)
{
	int nz;
	int n;
	int ic;
	complex s1;

	if(nin <= 2)
		n = nin;
	else
		n = 2;

	ic = 0;
	nz = 0;

	if(n >= 1)
	{
		s1 = y[0];
		nz = 1;
		y[0].re = 0.0;
		y[0].im = 0.0;

		if(-zr.re + log(rt_hypotd_snf(s1.re, s1.im)) >= -elim)
		{
			s1 = calccs(zr, s1, tol);

			if(cuchk(s1, ascle, tol) == 0)
			{
				y[0] = s1;
				nz = 0;
				ic = 1;
			}
		}
	}

	if(n == 2)
	{
		if(nz == 1)
			nz = 2;
		else
			nz = 1;

		s1 = y[1];
		y[1].re = 0.0;
		y[1].im = 0.0;

		if(-zr.re + log(rt_hypotd_snf(s1.re, s1.im)) >= -elim)
		{
			s1 = calccs(zr, s1, tol);

			if(cuchk(s1, ascle, tol) == 0)
			{
				y[1] = s1;
				nz = 0;
				ic = 2;
			}
		}

		if(ic < 2)
		{
			y[0].re = 0.0;
			y[0].im = 0.0;
			nz = 2;
		}
	}

	return nz;
}

complex calccs(complex zr, complex s1, double tol)
{
	complex cs;
	double tmp_re;
	double tmp_im;
	tmp_re = s1.re;
	tmp_im = s1.im;

	if ((s1.im == 0.0) && rtIsNaN(s1.re))
		;
	else if((fabs(s1.re) > 8.9884656743115785E+307) || (fabs(s1.im) > 8.9884656743115785E+307))
	{
		tmp_re = log(rt_hypotd_snf(s1.re / 2.0, s1.im / 2.0)) + 0.69314718055994529;
		tmp_im = rt_atan2d_snf(s1.im, s1.re);
	}
	else
	{
		tmp_re = log(rt_hypotd_snf(s1.re, s1.im));
		tmp_im = rt_atan2d_snf(s1.im, s1.re);
	}

	cs.im = -zr.im + tmp_im;
	tmp_re += -zr.re;
	cs.re = exp(tmp_re) / tol * cos(cs.im);
	cs.im = exp(tmp_re) / tol * sin(cs.im);
	return cs;
}

int ckscl(complex zr, complex *y, double elim)
{
	int nz;
	complex s1;
	s1 = *y;
	nz = 1;
	y->re = 0.0;
	y->im = 0.0;

	if(-zr.re + log(rt_hypotd_snf(s1.re, s1.im)) >= -elim)
	{
		s1 = calccs(zr, s1, 2.2204460492503131E-16);

		if(cuchk(s1, 1.0020841800044864E-289, 2.2204460492503131E-16) == 0)
		{
			*y = s1;
			nz = 0;
		}
	}

	return nz;
}

int b_cbknu(complex z, double fnu, int kode, int nin, complex y[2], double tol, double elim, double alim)
{
	int nz;
	int n;
	double yy;
	double caz;
	double cssr[3];
	double csrr[3];
	double bry[3];
	int iflag;
	double rz_re;
	double fks;
	double rz_im;
	double etest;
	int inu;
	double a1;
	double dnu;
	double rk;
	bool goto_mw110;
	double dnu2;
	double ak;
	double fhs;
	double s1_re;
	double s1_im;
	double s2_re;
	double s2_im;
	complex zd;
	double ck_re;
	double ck_im;
	int inub;
	bool goto_mw225;
	bool goto_mw240;
	bool goto_mw270;
	bool guard1 = false;
	bool guard2 = false;
	bool goto_mw210;
	complex cs;
	complex p2;
	double coef_re;
	complex fmu;
	double coef_im;
	int kflag;
	double a2;
	double t1;
	int i;
	bool exitg3;
	complex f;
	static const double dv10[8] = { 0.57721566490153287, -0.042002635034095237,
		-0.042197734555544333, 0.0072189432466631, -0.00021524167411495098,
		-2.0134854780788239E-5, 1.1330272319816959E-6, 6.1160951044814161E-9 };

	complex p1;
	int k;
	bool exitg2;
	long long i2;
	int j;
	complex cy[2];
	bool exitg1;
	double b_etest;
	bool b_guard1 = false;
	bool b_guard2 = false;

	if(nin <= 2)
		n = nin;
	else
		n = 2;

	yy = z.im;
	caz = rt_hypotd_snf(z.re, z.im);
	cssr[0] = 1.0 / tol;
	cssr[1] = 1.0;
	cssr[2] = tol;
	csrr[0] = tol;
	csrr[1] = 1.0;
	csrr[2] = 1.0 / tol;
	bry[0] = 2.2250738585072014E-305 / tol;
	bry[1] = tol / 2.2250738585072014E-305;
	bry[2] = 1.7976931348623157E+308;
	nz = 0;
	iflag = 0;

	if(z.im == 0.0)
	{
		rz_re = 2.0 / z.re;
		rz_im = 0.0;
	}
	else if (z.re == 0.0)
	{
		rz_re = 0.0;
		rz_im = -(2.0 / z.im);
	}
	else
	{
		fks = fabs(z.re);
		etest = fabs(z.im);

		if(fks > etest)
		{
			a1 = z.im / z.re;
			rk = z.re + a1 * z.im;
			rz_re = (2.0 + a1 * 0.0) / rk;
			rz_im = (0.0 - a1 * 2.0) / rk;
		}
		else if(etest == fks)
		{
			if(z.re > 0.0)
				rk = 0.5;
			else
				rk = -0.5;

			if(z.im > 0.0)
				a1 = 0.5;
			else
				a1 = -0.5;

			rz_re = (2.0 * rk + 0.0 * a1) / fks;
			rz_im = (0.0 * rk - 2.0 * a1) / fks;
		}
		else
		{
			a1 = z.re / z.im;
			rk = z.im + a1 * z.re;
			rz_re = a1 * 2.0 / rk;
			rz_im = (a1 * 0.0 - 2.0) / rk;
		}
	}

	inu = (int)(fnu + 0.5);
	dnu = fnu - (double)inu;
	goto_mw110 = (fabs(dnu) == 0.5);

	if((!goto_mw110) && (fabs(dnu) > tol))
		dnu2 = dnu * dnu;
	else
		dnu2 = 0.0;

	ak = 1.0;
	fhs = 0.0;
	s1_re = 0.0;
	s1_im = 0.0;
	s2_re = 0.0;
	s2_im = 0.0;
	zd.re = 0.0;
	zd.im = 0.0;
	ck_re = 1.0;
	ck_im = 0.0;
	inub = 1;
	goto_mw225 = false;
	goto_mw240 = false;
	goto_mw270 = false;

	if(goto_mw110 || (caz > 2.0))
		goto_mw110 = true;
	else
		goto_mw110 = false;

	guard1 = false;
	guard2 = false;

	if(!goto_mw110)
	{
		rk = 1.0;
		p2.re = rz_re;
		p2.im = rz_im;
		b_log(&p2);
		fmu.re = dnu * p2.re;
		fmu.im = dnu * p2.im;

		if(dnu != 0.0)
		{
			rk = dnu * 3.1415926535897931;
			rk /= sin(rk);
			cs = fmu;
			b_sinh(&cs);
			p2.re = 1.0 / dnu * cs.re;
			p2.im = 1.0 / dnu * cs.im;
		}

		etest = 1.0 + dnu;
		gammaln(&etest);
		a2 = exp(-etest);
		t1 = 1.0 / (a2 * rk);

		if(fabs(dnu) > 0.1)
			etest = (t1 - a2) / (dnu + dnu);
		else
		{
			a1 = 0.57721566490153287;
			i = 2;
			exitg3 = false;

			while((!exitg3) && (i < 9))
			{
				ak *= dnu2;
				s1_re = dv10[i - 1] * ak;
				a1 += s1_re;

				if(fabs(s1_re) < tol)
					exitg3 = true;
				else
					i++;
			}

			etest = -a1;
		}

		etest *= rk;
		cs = fmu;
		b_cosh(&cs);
		f.re = p2.re * (0.5 * (t1 + a2) * rk) + etest * cs.re;
		f.im = p2.im * (0.5 * (t1 + a2) * rk) + etest * cs.im;
		b_exp(&fmu);
		p1.re = 0.5 / a2 * fmu.re;
		p1.im = 0.5 / a2 * fmu.im;
		ak = 0.5 / t1;

		if(fmu.im == 0.0)
		{
			cs.re = ak / fmu.re;
			cs.im = 0.0;
		}
		else if (fmu.re == 0.0)
		{
			if(ak == 0.0)
			{
				cs.re = 0.0 / fmu.im;
				cs.im = 0.0;
			}
			else
			{
				cs.re = 0.0;
				cs.im = -(ak / fmu.im);
			}
		}
		else
		{
			fks = fabs(fmu.re);
			etest = fabs(fmu.im);

			if(fks > etest)
			{
				a1 = fmu.im / fmu.re;
				rk = fmu.re + a1 * fmu.im;
				cs.re = (ak + a1 * 0.0) / rk;
				cs.im = (0.0 - a1 * ak) / rk;
			}
			else if(etest == fks)
			{
				if(fmu.re > 0.0)
					rk = 0.5;
				else
					rk = -0.5;

				if(fmu.im > 0.0)
					a1 = 0.5;
				else
					a1 = -0.5;

				cs.re = (ak * rk + 0.0 * a1) / fks;
				cs.im = (0.0 * rk - ak * a1) / fks;
			}
			else
			{
				a1 = fmu.re / fmu.im;
				rk = fmu.im + a1 * fmu.re;
				cs.re = a1 * ak / rk;
				cs.im = (a1 * 0.0 - ak) / rk;
			}
		}

		s1_re = f.re;
		s1_im = f.im;
		s2_re = p1.re;
		s2_im = p1.im;
		ak = 1.0;
		a1 = 1.0;
		a2 = 1.0 - dnu2;

		if((inu > 0) || (n > 1))
			goto_mw110 = true;
		else
			goto_mw110 = false;

		if(!goto_mw110)
		{
			if(caz >= tol)
			{
				coef_re = 0.25 * (z.re * z.re - z.im * z.im);
				coef_im = 0.25 * (z.re * z.im + z.im * z.re);
				t1 = 0.25 * caz * caz;

				do
				{
					f.re *= ak;
					f.im *= ak;
					f.re += p1.re;
					f.im += p1.im;
					f.re += cs.re;
					f.im += cs.im;
					f.re *= 1.0 / a2;
					f.im *= 1.0 / a2;
					p1.re *= 1.0 / (ak - dnu);
					p1.im *= 1.0 / (ak - dnu);
					cs.re *= 1.0 / (ak + dnu);
					cs.im *= 1.0 / (ak + dnu);
					etest = ck_re;
					ck_re = 1.0 / ak * (ck_re * coef_re - ck_im * coef_im);
					ck_im = 1.0 / ak * (etest * coef_im + ck_im * coef_re);
					s1_re += ck_re * f.re - ck_im * f.im;
					s1_im += ck_re * f.im + ck_im * f.re;
					a1 = a1 * t1 / ak;
					a2 = ((a2 + ak) + ak) + 1.0;
					ak++;
				}
				while(!!(a1 > tol));
			}

			y[0].re = s1_re;
			y[0].im = s1_im;

			if(kode != 1)
			{
				y[0] = z;
				b_exp(&y[0]);
				etest = y[0].re;
				rk = y[0].im;
				y[0].re = etest * s1_re - rk * s1_im;
				y[0].im = etest * s1_im + rk * s1_re;
			}
		}
		else
		{
			if(caz >= tol)
			{
				coef_re = 0.25 * (z.re * z.re - z.im * z.im);
				coef_im = 0.25 * (z.re * z.im + z.im * z.re);
				t1 = 0.25 * caz * caz;

				do
				{
					f.re *= ak;
					f.im *= ak;
					f.re += p1.re;
					f.im += p1.im;
					f.re += cs.re;
					f.im += cs.im;
					f.re *= 1.0 / a2;
					f.im *= 1.0 / a2;
					p1.re *= 1.0 / (ak - dnu);
					p1.im *= 1.0 / (ak - dnu);
					cs.re *= 1.0 / (ak + dnu);
					cs.im *= 1.0 / (ak + dnu);
					etest = ck_re;
					ck_re = 1.0 / ak * (ck_re * coef_re - ck_im * coef_im);
					ck_im = 1.0 / ak * (etest * coef_im + ck_im * coef_re);
					s1_re += ck_re * f.re - ck_im * f.im;
					s1_im += ck_re * f.im + ck_im * f.re;
					etest = p1.re - f.re * ak;
					rk = p1.im - f.im * ak;
					s2_re += etest * ck_re - rk * ck_im;
					s2_im += etest * ck_im + rk * ck_re;
					a1 = a1 * t1 / ak;
					a2 = ((a2 + ak) + ak) + 1.0;
					ak++;
				}
				while(!!(a1 > tol));
			}

			kflag = 1;

			if((fnu + 1.0) * fabs(p2.re) > alim)
				kflag = 2;

			etest = cssr[kflag] * s2_re;
			s2_im *= cssr[kflag];
			s2_re = etest * rz_re - s2_im * rz_im;
			s2_im = etest * rz_im + s2_im * rz_re;
			s1_re *= cssr[kflag];
			s1_im *= cssr[kflag];

			if(kode != 1)
			{
				f = z;
				b_exp(&f);
				etest = s1_re;
				s1_re = s1_re * f.re - s1_im * f.im;
				s1_im = etest * f.im + s1_im * f.re;
				etest = s2_re;
				s2_re = s2_re * f.re - s2_im * f.im;
				s2_im = etest * f.im + s2_im * f.re;
			}

			goto_mw210 = true;
			guard1 = true;
		}
	}
	else
	{
		goto_mw210 = false;
		cs = z;
		b_sqrt(&cs);

		if(cs.im == 0.0)
		{
			coef_re = 1.2533141373155001 / cs.re;
			coef_im = 0.0;
		}
		else if(cs.re == 0.0)
		{
			coef_re = 0.0;
			coef_im = -(1.2533141373155001 / cs.im);
		}
		else
		{
			fks = fabs(cs.re);
			etest = fabs(cs.im);

			if(fks > etest)
			{
				a1 = cs.im / cs.re;
				rk = cs.re + a1 * cs.im;
				coef_re = (1.2533141373155001 + a1 * 0.0) / rk;
				coef_im = (0.0 - a1 * 1.2533141373155001) / rk;
			}
			else if(etest == fks)
			{
				if(cs.re > 0.0)
					rk = 0.5;
				else
					rk = -0.5;

				if(cs.im > 0.0)
					a1 = 0.5;
				else
					a1 = -0.5;

				coef_re = (1.2533141373155001 * rk + 0.0 * a1) / fks;
				coef_im = (0.0 * rk - 1.2533141373155001 * a1) / fks;
			}
			else
			{
				a1 = cs.re / cs.im;
				rk = cs.im + a1 * cs.re;
				coef_re = a1 * 1.2533141373155001 / rk;
				coef_im = (a1 * 0.0 - 1.2533141373155001) / rk;
			}
		}

		kflag = 1;

		if(!(kode == 2))
		{
			if(z.re > alim)
				iflag = 1;
			else
			{
				a1 = exp(-z.re);
				fmu.re = cos(z.im) * a1;
				fmu.im = sin(z.im) * -a1;
				a2 = coef_re;
				coef_re = coef_re * fmu.re - coef_im * fmu.im;
				coef_im = a2 * fmu.im + coef_im * fmu.re;
			}
		}

		if(fabs(dnu) == 0.5)
		{
			s1_re = coef_re;
			s1_im = coef_im;
			s2_re = coef_re;
			s2_im = coef_im;
			goto_mw210 = true;
		}
		else
		{
			ak = fabs(cos(3.1415926535897931 * dnu));

			if(ak == 0.0)
			{
				s1_re = coef_re;
				s1_im = coef_im;
				s2_re = coef_re;
				s2_im = coef_im;
				goto_mw210 = true;
			}
			else
			{
				fhs = fabs(0.25 - dnu2);

				if(fhs == 0.0)
				{
					s1_re = coef_re;
					s1_im = coef_im;
					s2_re = coef_re;
					s2_im = coef_im;
					goto_mw210 = true;
				}
			}
		}

		if(!goto_mw210)
		{
			if(z.re != 0.0)
				t1 = fabs(atan(z.im / z.re));
			else
				t1 = 1.5707963267948966;

			if(28.66666665740641 > caz)
			{
				ak = 1.8976999933151775 * ak / (tol * sqrt(sqrt(caz)));
				ak = (log(ak) + caz * cos(3.0 * t1 / (1.0 + caz)) / (1.0 + 0.008 * caz)) / cos(14.7 * t1 / (28.0 + caz));
				ak = 0.12125 * ak * ak / caz + 1.5;
				guard2 = true;
			}
			else
			{
				etest = ak / (3.1415926535897931 * caz * tol);
				ak = 1.0;

				if(etest >= 1.0)
				{
					fks = 2.0;
					rk = (caz + caz) + 2.0;
					a1 = 0.0;
					a2 = 1.0;
					goto_mw110 = true;
					i = 1;
					exitg2 = false;

					while((!exitg2) && (i < 31))
					{
						s1_re = a2;
						a2 = rk / (ak + 1.0) * a2 - fhs / fks * a1;
						a1 = s1_re;
						rk += 2.0;
						fks = ((fks + ak) + ak) + 2.0;
						fhs = (fhs + ak) + ak;
						ak++;

						if(etest < fabs(a2) * ak)
						{
							goto_mw110 = false;
							exitg2 = true;
						}
						else
							i++;
					}

					if(goto_mw110)
						nz = -2;
					else
					{
						ak += 1.909859317102744 * t1 * sqrt(28.66666665740641 / caz);
						fhs = fabs(0.25 - dnu2);
						guard2 = true;
					}
				}
				else
					guard2 = true;
			}
		}
		else
			guard1 = true;
	}

	if(guard2)
	{
		b_fix(&ak);
		k = (int)ak;
		ak = k;
		fks = (double)k * (double)k;
		p1.re = 0.0;
		p1.im = 0.0;
		p2.re = tol;
		p2.im = 0.0;
		cs = p2;

		for(i = 1; i <= k; i++)
		{
			a1 = fks - ak;
			rk = 2.0 / (ak + 1.0);
			fmu = p2;
			a2 = (ak + z.re) * rk;
			etest = yy * rk;
			rk = p2.re;
			p2.re = p2.re * a2 - p2.im * etest;
			p2.im = rk * etest + p2.im * a2;
			p2.re -= p1.re;
			p2.im -= p1.im;
			p2.re *= (fks + ak) / (a1 + fhs);
			p2.im *= (fks + ak) / (a1 + fhs);
			p1 = fmu;
			cs.re += p2.re;
			cs.im += p2.im;
			fks = (a1 - ak) + 1.0;
			ak--;
		}

		fmu.re = 1.0 / rt_hypotd_snf(cs.re, cs.im);
		cs.im = -cs.im;
		etest = cs.re;
		cs.re = cs.re * fmu.re - cs.im * 0.0;
		cs.im = etest * 0.0 + cs.im * fmu.re;
		rk = fmu.re * p2.re - 0.0 * p2.im;
		etest = fmu.re * p2.im + 0.0 * p2.re;
		a2 = coef_re * rk - coef_im * etest;
		coef_im = coef_re * etest + coef_im * rk;
		s1_re = a2 * cs.re - coef_im * cs.im;
		s1_im = a2 * cs.im + coef_im * cs.re;

		if((inu > 0) || (n > 1))
			goto_mw110 = true;
		else
			goto_mw110 = false;

		if(!goto_mw110)
		{
			zd = z;

			if(iflag == 1)
				goto_mw270 = true;
			else
				goto_mw240 = true;
		}
		else
		{
			fmu.re = 1.0 / rt_hypotd_snf(p2.re, p2.im);
			etest = p1.re;
			p1.re = p1.re * fmu.re - p1.im * 0.0;
			p1.im = etest * 0.0 + p1.im * fmu.re;
			p2.im = -p2.im;
			rk = p2.re;
			p2.re = p2.re * fmu.re - p2.im * 0.0;
			p2.im = rk * 0.0 + p2.im * fmu.re;
			rk = p1.re * p2.im + p1.im * p2.re;
			ak = (dnu + 0.5) - (p1.re * p2.re - p1.im * p2.im);
			a2 = 0.0 - (p1.re * p2.im + p1.im * p2.re);

			if(z.im == 0.0)
			{
				if(0.0 - rk == 0.0)
				{
					etest = ak / z.re;
					rk = 0.0;
				}
				else if(ak == 0.0)
				{
					etest = 0.0;
					rk = (0.0 - rk) / z.re;
				}
				else
				{
					etest = ak / z.re;
					rk = (0.0 - rk) / z.re;
				}
			}
			else if(z.re == 0.0)
			{
				if(ak == 0.0)
				{
					etest = (0.0 - rk) / z.im;
					rk = 0.0;
				}
				else if(0.0 - rk == 0.0)
				{
					etest = 0.0;
					rk = -(ak / z.im);
				}
				else
				{
					etest = (0.0 - rk) / z.im;
					rk = -(ak / z.im);
				}
			}
			else
			{
				fks = fabs(z.re);
				etest = fabs(z.im);

				if(fks > etest)
				{
					a1 = z.im / z.re;
					rk = z.re + a1 * z.im;
					etest = (ak + a1 * a2) / rk;
					rk = (a2 - a1 * ak) / rk;
				}
				else if(etest == fks)
				{
					if(z.re > 0.0)
						rk = 0.5;
					else
						rk = -0.5;

					if(z.im > 0.0)
						a1 = 0.5;
					else
						a1 = -0.5;

					etest = (ak * rk + a2 * a1) / fks;
					rk = (a2 * rk - ak * a1) / fks;
				}
				else
				{
					a1 = z.re / z.im;
					rk = z.im + a1 * z.re;
					etest = (a1 * ak + a2) / rk;
					rk = (a1 * a2 - ak) / rk;
				}
			}

			etest++;
			s2_re = s1_re * etest - s1_im * rk;
			s2_im = s1_re * rk + s1_im * etest;
			goto_mw210 = true;
		}

		guard1 = true;
	}

	if(guard1)
	{
		if(goto_mw240 || goto_mw270)
			;
		else if(goto_mw210)
		{
			ck_re = (dnu + 1.0) * rz_re;
			ck_im = (dnu + 1.0) * rz_im;

			if(n == 1)
				inu--;

			if(inu > 0)
			{
				if(iflag == 1)
				{
					a2 = 0.5 * elim;
					coef_re = exp(-elim);
					zd = z;
					rk = z.re;
					k = inu;
					j = 1;

					for(i = 0; i < 2; i++)
					{
						cy[i].re = 0.0;
						cy[i].im = 0.0;
					}

					i = 1;
					exitg1 = false;

					while((!exitg1) && (i <= inu))
					{
						f.re = s2_re;
						f.im = s2_im;
						etest = s2_re;
						s2_re = (s2_re * ck_re - s2_im * ck_im) + s1_re;
						s2_im = (etest * ck_im + s2_im * ck_re) + s1_im;
						s1_re = f.re;
						s1_im = f.im;
						ck_re += rz_re;
						ck_im += rz_im;
						etest = log(rt_hypotd_snf(s2_re, s2_im));
						b_guard1 = false;
						b_guard2 = false;

						if(-rk + etest >= -elim)
						{
							cs.re = s2_re;
							cs.im = s2_im;
							b_log(&cs);
							p2.re = cs.re + -zd.re;
							p2.im = cs.im + -zd.im;
							p1.re = exp(p2.re) / tol * cos(p2.im);
							p1.im = exp(p2.re) / tol * sin(p2.im);

							if(cuchk(p1, bry[0], tol) == 0)
							{
								j = 1 - j;
								cy[j] = p1;

								if(k == i - 1)
								{
									kflag = 0;
									inub = i + 1;
									s2_re = cy[j].re;
									s2_im = cy[j].im;
									j = -j;
									s1_re = cy[j + 1].re;
									s1_im = cy[j + 1].im;

									if(inub <= inu)
										goto_mw225 = true;
									else
									{
										if(n == 1)
										{
											s1_re = s2_re;
											s1_im = s2_im;
										}

										goto_mw240 = true;
									}

									exitg1 = true;
								}
								else
								{
									k = i;
									b_guard1 = true;
								}
							}
							else
								b_guard2 = true;
						}
						else
							b_guard2 = true;

						if(b_guard2)
						{
							if(etest >= a2)
							{
								rk -= elim;
								s1_re = f.re * coef_re - f.im * 0.0;
								s1_im = f.re * 0.0 + f.im * coef_re;
								etest = s2_re;
								s2_re = s2_re * coef_re - s2_im * 0.0;
								s2_im = etest * 0.0 + s2_im * coef_re;
								zd.re = rk;
								zd.im = yy;
							}

							b_guard1 = true;
						}

						if(b_guard1)
							i++;
					}

					if(goto_mw225 || goto_mw240)
						goto_mw110 = true;
					else
						goto_mw110 = false;

					if(!goto_mw110)
					{
						if(n == 1)
						{
							s1_re = s2_re;
							s1_im = s2_im;
						}

						goto_mw270 = true;
					}
				}
				else
					goto_mw225 = true;
			}

			if(goto_mw225 || goto_mw240 || goto_mw270)
				goto_mw110 = true;
			else
				goto_mw110 = false;

			if(!goto_mw110)
			{
				if(n == 1)
				{
					s1_re = s2_re;
					s1_im = s2_im;
				}

				zd = z;

				if(iflag == 1)
					;
				else
					goto_mw240 = true;
			}
		}
		else
			goto_mw225 = true;

		if(goto_mw225 || goto_mw240)
		{
			if(goto_mw225)
			{
				p1.re = csrr[kflag];
				a2 = bry[kflag];

				while(inub <= inu)
				{
					f.re = s2_re;
					f.im = s2_im;
					etest = s2_re;
					s2_re = (s2_re * ck_re - s2_im * ck_im) + s1_re;
					s2_im = (etest * ck_im + s2_im * ck_re) + s1_im;
					s1_re = f.re;
					s1_im = f.im;
					ck_re += rz_re;
					ck_im += rz_im;

					if(kflag + 1 < 3)
					{
						p2.re = s2_re * p1.re - s2_im * 0.0;
						p2.im = s2_re * 0.0 + s2_im * p1.re;
						etest = fabs(p2.re);
						rk = fabs(p2.im);

						if((etest >= rk) || rtIsNaN(rk))
							b_etest = etest;
						else
							b_etest = rk;

						if(b_etest > a2)
						{
							kflag++;
							a2 = bry[kflag];
							s1_re = cssr[kflag] * (f.re * p1.re - f.im * 0.0);
							s1_im = cssr[kflag] * (f.re * 0.0 + f.im * p1.re);
							s2_re = cssr[kflag] * p2.re;
							s2_im = cssr[kflag] * p2.im;
							p1.re = csrr[kflag];
						}
					}

					inub++;
				}

				if(n == 1)
				{
					s1_re = s2_re;
					s1_im = s2_im;
				}
			}

			y[0].re = csrr[kflag] * s1_re;
			y[0].im = csrr[kflag] * s1_im;

			if(n == 2)
			{
				y[1].re = csrr[kflag] * s2_re;
				y[1].im = csrr[kflag] * s2_im;
			}
		}
		else
		{
			y[0].re = s1_re;
			y[0].im = s1_im;

			if(n != 1)
			{
				y[1].re = s2_re;
				y[1].im = s2_im;
			}

			nz = b_ckscl(zd, n, y, bry[0], tol, elim);

			if(n <= nz)
				;
			else
			{
				y[nz].re *= tol;
				y[nz].im *= tol;
				i2 = (long long)n - nz;

				if(i2 > 2147483647LL)
					i2 = 2147483647LL;
				else
				{
					if(i2 < -2147483648LL)
						i2 = -2147483648LL;
				}

				if((int)i2 >= 2)
				{
					k = nz + 1;
					y[k].re *= tol;
					y[k].im *= tol;
				}
			}
		}
	}

	return nz;
}

void cbknu(complex z, double fnu, int kode, double alim, complex *y, int *nz)
{
	double yy;
	double caz;
	int iflag;
	double rz_re;
	double fks;
	double rz_im;
	double etest;
	int inu;
	double dnu;
	bool goto_mw110;
	double dnu2;
	double tm;
	double ak;
	double fhs;
	complex s1;
	double s2_re;
	double s2_im;
	complex zd;
	double ck_re;
	double ck_im;
	int inub;
	bool goto_mw225;
	bool goto_mw240;
	bool goto_mw270;
	bool guard1 = false;
	bool guard2 = false;
	double a2;
	complex coef;
	complex smu;
	complex fmu;
	int kflag;
	double rk;
	double t1;
	int i;
	bool exitg3;
	complex p1;
	static const double dv6[8] = { 0.57721566490153287, -0.042002635034095237,
		-0.042197734555544333, 0.0072189432466631, -0.00021524167411495098,
		-2.0134854780788239E-5, 1.1330272319816959E-6, 6.1160951044814161E-9 };

	double p2_re;
	int k;
	double p2_im;
	double cs_re;
	double a1;
	static const double dv7[3] = { 2.2204460492503131E-16, 1.0, 4.503599627370496E+15 };

	static const double dv8[3] = { 1.0020841800044864E-289, 9.9792015476736E+288, 1.7976931348623157E+308 };

	bool earlyExit;
	bool exitg2;
	long long i1;
	int j;
	complex cy[2];
	bool exitg1;
	static const double dv9[3] = { 4.503599627370496E+15, 1.0, 2.2204460492503131E-16 };

	double b_etest;
	bool b_guard1 = false;
	bool b_guard2 = false;
	y->re = 0.0;
	y->im = 0.0;
	yy = z.im;
	caz = rt_hypotd_snf(z.re, z.im);
	*nz = 0;
	iflag = 0;

	if(z.im == 0.0)
	{
		rz_re = 2.0 / z.re;
		rz_im = 0.0;
	}
	else if(z.re == 0.0)
	{
		rz_re = 0.0;
		rz_im = -(2.0 / z.im);
	}
	else
	{
		fks = fabs(z.re);
		etest = fabs(z.im);

		if(fks > etest)
		{
			fks = z.im / z.re;
			etest = z.re + fks * z.im;
			rz_re = (2.0 + fks * 0.0) / etest;
			rz_im = (0.0 - fks * 2.0) / etest;
		}
		else if(etest == fks)
		{
			if(z.re > 0.0)
				etest = 0.5;
			else
				etest = -0.5;

			if(z.im > 0.0)
				tm = 0.5;
			else
				tm = -0.5;

			rz_re = (2.0 * etest + 0.0 * tm) / fks;
			rz_im = (0.0 * etest - 2.0 * tm) / fks;
		}
		else
		{
			fks = z.re / z.im;
			etest = z.im + fks * z.re;
			rz_re = fks * 2.0 / etest;
			rz_im = (fks * 0.0 - 2.0) / etest;
		}
	}

	inu = (int)(fnu + 0.5);
	dnu = fnu - (double)inu;
	goto_mw110 = (fabs(dnu) == 0.5);

	if((!goto_mw110) && (fabs(dnu) > 2.2204460492503131E-16))
		dnu2 = dnu * dnu;
	else
		dnu2 = 0.0;

	ak = 1.0;
	fhs = 0.0;
	s1.re = 0.0;
	s1.im = 0.0;
	s2_re = 0.0;
	s2_im = 0.0;
	zd.re = 0.0;
	zd.im = 0.0;
	ck_re = 1.0;
	ck_im = 0.0;
	inub = 1;
	goto_mw225 = false;
	goto_mw240 = false;
	goto_mw270 = false;

	if(goto_mw110 || (caz > 2.0))
		goto_mw110 = true;
	else
		goto_mw110 = false;

	guard1 = false;
	guard2 = false;

	if(!goto_mw110)
	{
		a2 = 1.0;
		smu.re = rz_re;
		smu.im = rz_im;
		b_log(&smu);
		fmu.re = dnu * smu.re;
		fmu.im = dnu * smu.im;

		if(dnu != 0.0)
		{
			a2 = dnu * 3.1415926535897931;
			a2 /= sin(a2);
			coef = fmu;
			b_sinh(&coef);
			smu.re = 1.0 / dnu * coef.re;
			smu.im = 1.0 / dnu * coef.im;
		}

		etest = 1.0 + dnu;
		gammaln(&etest);
		rk = exp(-etest);
		t1 = 1.0 / (rk * a2);

		if(fabs(dnu) > 0.1)
			etest = (t1 - rk) / (dnu + dnu);
		else
		{
			fks = 0.57721566490153287;
			i = 2;
			exitg3 = false;

			while((!exitg3) && (i < 9))
			{
				ak *= dnu2;
				tm = dv6[i - 1] * ak;
				fks += tm;

				if(fabs(tm) < 2.2204460492503131E-16)
					exitg3 = true;
				else
					i++;
			}

			etest = -fks;
		}

		etest *= a2;
		coef = fmu;
		b_cosh(&coef);
		p1.re = smu.re * (0.5 * (t1 + rk) * a2) + etest * coef.re;
		p1.im = smu.im * (0.5 * (t1 + rk) * a2) + etest * coef.im;
		b_exp(&fmu);
		p2_re = 0.5 / rk * fmu.re;
		p2_im = 0.5 / rk * fmu.im;
		ak = 0.5 / t1;

		if(fmu.im == 0.0)
		{
			cs_re = ak / fmu.re;
			tm = 0.0;
		}
		else if(fmu.re == 0.0)
		{
			if(ak == 0.0)
			{
				cs_re = 0.0 / fmu.im;
				tm = 0.0;
			}
			else
			{
				cs_re = 0.0;
				tm = -(ak / fmu.im);
			}
		}
		else
		{
			fks = fabs(fmu.re);
			etest = fabs(fmu.im);

			if(fks > etest)
			{
				fks = fmu.im / fmu.re;
				etest = fmu.re + fks * fmu.im;
				cs_re = (ak + fks * 0.0) / etest;
				tm = (0.0 - fks * ak) / etest;
			}
			else if(etest == fks)
			{
				if(fmu.re > 0.0)
					etest = 0.5;
				else
					etest = -0.5;

				if(fmu.im > 0.0)
					tm = 0.5;
				else
					tm = -0.5;

				cs_re = (ak * etest + 0.0 * tm) / fks;
				tm = (0.0 * etest - ak * tm) / fks;
			}
			else
			{
				fks = fmu.re / fmu.im;
				etest = fmu.im + fks * fmu.re;
				cs_re = fks * ak / etest;
				tm = (fks * 0.0 - ak) / etest;
			}
		}

		s1 = p1;
		s2_re = p2_re;
		s2_im = p2_im;
		ak = 1.0;
		a1 = 1.0;
		rk = 1.0 - dnu2;

		if(!(inu > 0))
		{
			if(caz >= 2.2204460492503131E-16)
			{
				coef.re = 0.25 * (z.re * z.re - z.im * z.im);
				coef.im = 0.25 * (z.re * z.im + z.im * z.re);
				t1 = 0.25 * caz * caz;

				do
				{
					p1.re *= ak;
					p1.im *= ak;
					p1.re += p2_re;
					p1.im += p2_im;
					p1.re += cs_re;
					p1.im += tm;
					p1.re *= 1.0 / rk;
					p1.im *= 1.0 / rk;
					p2_re *= 1.0 / (ak - dnu);
					p2_im *= 1.0 / (ak - dnu);
					cs_re *= 1.0 / (ak + dnu);
					tm *= 1.0 / (ak + dnu);
					etest = ck_re;
					ck_re = 1.0 / ak * (ck_re * coef.re - ck_im * coef.im);
					ck_im = 1.0 / ak * (etest * coef.im + ck_im * coef.re);
					s1.re += ck_re * p1.re - ck_im * p1.im;
					s1.im += ck_re * p1.im + ck_im * p1.re;
					a1 = a1 * t1 / ak;
					rk = ((rk + ak) + ak) + 1.0;
					ak++;
				}
				while(!!(a1 > 2.2204460492503131E-16));
			}

			*y = s1;

			if(kode != 1)
			{
				coef = z;
				b_exp(&coef);
				y->re = coef.re * s1.re - coef.im * s1.im;
				y->im = coef.re * s1.im + coef.im * s1.re;
			}
		}
		else
		{
			if(caz >= 2.2204460492503131E-16)
			{
				coef.re = 0.25 * (z.re * z.re - z.im * z.im);
				coef.im = 0.25 * (z.re * z.im + z.im * z.re);
				t1 = 0.25 * caz * caz;

				do
				{
					p1.re *= ak;
					p1.im *= ak;
					p1.re += p2_re;
					p1.im += p2_im;
					p1.re += cs_re;
					p1.im += tm;
					p1.re *= 1.0 / rk;
					p1.im *= 1.0 / rk;
					p2_re *= 1.0 / (ak - dnu);
					p2_im *= 1.0 / (ak - dnu);
					cs_re *= 1.0 / (ak + dnu);
					tm *= 1.0 / (ak + dnu);
					etest = ck_re;
					ck_re = 1.0 / ak * (ck_re * coef.re - ck_im * coef.im);
					ck_im = 1.0 / ak * (etest * coef.im + ck_im * coef.re);
					s1.re += ck_re * p1.re - ck_im * p1.im;
					s1.im += ck_re * p1.im + ck_im * p1.re;
					a2 = p2_re - p1.re * ak;
					etest = p2_im - p1.im * ak;
					s2_re += a2 * ck_re - etest * ck_im;
					s2_im += a2 * ck_im + etest * ck_re;
					a1 = a1 * t1 / ak;
					rk = ((rk + ak) + ak) + 1.0;
					ak++;
				}
				while(!!(a1 > 2.2204460492503131E-16));
			}

			kflag = 1;

			if((fnu + 1.0) * fabs(smu.re) > alim)
				kflag = 2;

			etest = dv9[kflag] * s2_re;
			s2_im *= dv9[kflag];
			s2_re = etest * rz_re - s2_im * rz_im;
			s2_im = etest * rz_im + s2_im * rz_re;
			s1.re *= dv9[kflag];
			s1.im *= dv9[kflag];

			if(kode != 1)
			{
				p1 = z;
				b_exp(&p1);
				etest = s1.re;
				s1.re = s1.re * p1.re - s1.im * p1.im;
				s1.im = etest * p1.im + s1.im * p1.re;
				etest = s2_re;
				s2_re = s2_re * p1.re - s2_im * p1.im;
				s2_im = etest * p1.im + s2_im * p1.re;
			}

			goto_mw110 = true;
			guard1 = true;
		}
	}
	else
	{
		goto_mw110 = false;
		coef = z;
		b_sqrt(&coef);

		if(coef.im == 0.0)
		{
			coef.re = 1.2533141373155001 / coef.re;
			coef.im = 0.0;
		}
		else if(coef.re == 0.0)
		{
			coef.re = 0.0;
			coef.im = -(1.2533141373155001 / coef.im);
		}
		else
		{
			fks = fabs(coef.re);
			etest = fabs(coef.im);

			if(fks > etest)
			{
				fks = coef.im / coef.re;
				etest = coef.re + fks * coef.im;
				coef.re = (1.2533141373155001 + fks * 0.0) / etest;
				coef.im = (0.0 - fks * 1.2533141373155001) / etest;
			}
			else if(etest == fks)
			{
				if(coef.re > 0.0)
					etest = 0.5;
				else
					etest = -0.5;

				if(coef.im > 0.0)
					tm = 0.5;
				else
					tm = -0.5;

				coef.re = (1.2533141373155001 * etest + 0.0 * tm) / fks;
				coef.im = (0.0 * etest - 1.2533141373155001 * tm) / fks;
			}
			else
			{
				fks = coef.re / coef.im;
				etest = coef.im + fks * coef.re;
				coef.re = fks * 1.2533141373155001 / etest;
				coef.im = (fks * 0.0 - 1.2533141373155001) / etest;
			}
		}

		kflag = 1;

		if(!(kode == 2))
		{
			if(z.re > alim)
				iflag = 1;
			else
			{
				etest = exp(-z.re);
				fmu.re = cos(z.im) * etest;
				fmu.im = sin(z.im) * -etest;
				rk = coef.re;
				coef.re = coef.re * fmu.re - coef.im * fmu.im;
				coef.im = rk * fmu.im + coef.im * fmu.re;
			}
		}

		if(fabs(dnu) == 0.5)
		{
			s1 = coef;
			s2_re = coef.re;
			s2_im = coef.im;
			goto_mw110 = true;
		}
		else
		{
			ak = fabs(cos(3.1415926535897931 * dnu));

			if(ak == 0.0)
			{
				s1 = coef;
				s2_re = coef.re;
				s2_im = coef.im;
				goto_mw110 = true;
			}
			else
			{
				fhs = fabs(0.25 - dnu2);

				if(fhs == 0.0)
				{
					s1 = coef;
					s2_re = coef.re;
					s2_im = coef.im;
					goto_mw110 = true;
				}
			}
		}

		if(!goto_mw110)
		{
			if(z.re != 0.0)
				t1 = fabs(atan(z.im / z.re));
			else
				t1 = 1.5707963267948966;

			if(28.66666665740641 > caz)
			{
				ak = 1.8976999933151775 * ak / (2.2204460492503131E-16 * sqrt(sqrt(caz)));
				ak = (log(ak) + caz * cos(3.0 * t1 / (1.0 + caz)) / (1.0 + 0.008 * caz)) / cos(14.7 * t1 / (28.0 + caz));
				ak = 0.12125 * ak * ak / caz + 1.5;
				guard2 = true;
			}
			else
			{
				etest = ak / (3.1415926535897931 * caz * 2.2204460492503131E-16);
				ak = 1.0;

				if(etest >= 1.0)
				{
					fks = 2.0;
					rk = (caz + caz) + 2.0;
					a1 = 0.0;
					a2 = 1.0;
					earlyExit = true;
					i = 1;
					exitg2 = false;

					while((!exitg2) && (i < 31))
					{
						tm = a2;
						a2 = rk / (ak + 1.0) * a2 - fhs / fks * a1;
						a1 = tm;
						rk += 2.0;
						fks = ((fks + ak) + ak) + 2.0;
						fhs = (fhs + ak) + ak;
						ak++;

						if(etest < fabs(a2) * ak)
						{
							earlyExit = false;
							exitg2 = true;
						}
						else
							i++;
					}

					if(earlyExit)
						*nz = -2;
					else
					{
						ak += 1.909859317102744 * t1 * sqrt(28.66666665740641 / caz);
						fhs = fabs(0.25 - dnu2);
						guard2 = true;
					}
				}
				else
					guard2 = true;
			}
		}
		else
			guard1 = true;
	}

	if(guard2)
	{
		b_fix(&ak);
		k = (int)ak;
		ak = k;
		fks = (double)k * (double)k;
		p1.re = 0.0;
		p1.im = 0.0;
		p2_re = 2.2204460492503131E-16;
		p2_im = 0.0;
		cs_re = 2.2204460492503131E-16;
		tm = 0.0;

		for(i = 1; i <= k; i++)
		{
			a1 = fks - ak;
			rk = 2.0 / (ak + 1.0);
			fmu.re = p2_re;
			fmu.im = p2_im;
			etest = (ak + z.re) * rk;
			rk *= yy;
			a2 = p2_re;
			p2_re = (fks + ak) / (a1 + fhs) * ((p2_re * etest - p2_im * rk) - p1.re);
			p2_im = (fks + ak) / (a1 + fhs) * ((a2 * rk + p2_im * etest) - p1.im);
			p1 = fmu;
			cs_re += p2_re;
			tm += p2_im;
			fks = (a1 - ak) + 1.0;
			ak--;
		}

		fmu.re = 1.0 / rt_hypotd_snf(cs_re, tm);
		tm = -tm;
		etest = cs_re;
		cs_re = cs_re * fmu.re - tm * 0.0;
		tm = etest * 0.0 + tm * fmu.re;
		etest = fmu.re * p2_re - 0.0 * p2_im;
		a2 = fmu.re * p2_im + 0.0 * p2_re;
		rk = coef.re * etest - coef.im * a2;
		etest = coef.re * a2 + coef.im * etest;
		s1.re = rk * cs_re - etest * tm;
		s1.im = rk * tm + etest * cs_re;

		if(!(inu > 0))
		{
			zd = z;

			if(iflag == 1)
				goto_mw270 = true;
			else
				goto_mw240 = true;
		}
		else
		{
			fmu.re = 1.0 / rt_hypotd_snf(p2_re, p2_im);
			etest = p1.re;
			p1.re = p1.re * fmu.re - p1.im * 0.0;
			p1.im = etest * 0.0 + p1.im * fmu.re;
			p2_im = -p2_im;
			a2 = p2_re;
			p2_re = p2_re * fmu.re - p2_im * 0.0;
			p2_im = a2 * 0.0 + p2_im * fmu.re;
			etest = p1.re * p2_im + p1.im * p2_re;
			ak = (dnu + 0.5) - (p1.re * p2_re - p1.im * p2_im);
			a2 = 0.0 - (p1.re * p2_im + p1.im * p2_re);

			if(z.im == 0.0)
			{
				if(0.0 - etest == 0.0)
				{
					rk = ak / z.re;
					etest = 0.0;
				}
				else if(ak == 0.0)
				{
					rk = 0.0;
					etest = (0.0 - etest) / z.re;
				}
				else
				{
					rk = ak / z.re;
					etest = (0.0 - etest) / z.re;
				}
			}
			else if(z.re == 0.0)
			{
				if(ak == 0.0)
				{
					rk = (0.0 - etest) / z.im;
					etest = 0.0;
				}
				else if(0.0 - etest == 0.0)
				{
					rk = 0.0;
					etest = -(ak / z.im);
				}
				else
				{
					rk = (0.0 - etest) / z.im;
					etest = -(ak / z.im);
				}
			}
			else
			{
				fks = fabs(z.re);
				etest = fabs(z.im);

				if(fks > etest)
				{
					fks = z.im / z.re;
					etest = z.re + fks * z.im;
					rk = (ak + fks * a2) / etest;
					etest = (a2 - fks * ak) / etest;
				}
				else if(etest == fks)
				{
					if(z.re > 0.0)
						etest = 0.5;
					else
						etest = -0.5;

					if(z.im > 0.0)
						tm = 0.5;
					else
						tm = -0.5;

					rk = (ak * etest + a2 * tm) / fks;
					etest = (a2 * etest - ak * tm) / fks;
				}
				else
				{
					fks = z.re / z.im;
					etest = z.im + fks * z.re;
					rk = (fks * ak + a2) / etest;
					etest = (fks * a2 - ak) / etest;
				}
			}

			rk++;
			s2_re = s1.re * rk - s1.im * etest;
			s2_im = s1.re * etest + s1.im * rk;
			goto_mw110 = true;
		}

		guard1 = true;
	}

	if(guard1)
	{
		if(goto_mw240 || goto_mw270)
			;
		else if(goto_mw110)
		{
			ck_re = (dnu + 1.0) * rz_re;
			ck_im = (dnu + 1.0) * rz_im;
			inu--;

			if(inu > 0)
			{
				if(iflag == 1)
				{
					zd = z;
					a2 = z.re;
					k = inu;
					j = 1;

					for(i = 0; i < 2; i++)
					{
						cy[i].re = 0.0;
						cy[i].im = 0.0;
					}

					i = 1;
					exitg1 = false;

					while((!exitg1) && (i <= inu))
					{
						cs_re = s2_re;
						tm = s2_im;
						etest = s2_re;
						s2_re = (s2_re * ck_re - s2_im * ck_im) + s1.re;
						s2_im = (etest * ck_im + s2_im * ck_re) + s1.im;
						s1.re = cs_re;
						s1.im = tm;
						ck_re += rz_re;
						ck_im += rz_im;
						etest = log(rt_hypotd_snf(s2_re, s2_im));
						b_guard1 = false;
						b_guard2 = false;

						if(-a2 + etest >= -700.92179369444591)
						{
							coef.re = s2_re;
							coef.im = s2_im;
							b_log(&coef);
							p2_re = coef.re + -zd.re;
							p2_im = coef.im + -zd.im;
							p1.re = exp(p2_re) / 2.2204460492503131E-16 * cos(p2_im);
							p1.im = exp(p2_re) / 2.2204460492503131E-16 * sin(p2_im);

							if(cuchk(p1, 1.0020841800044864E-289, 2.2204460492503131E-16) == 0)
							{
								j = 1 - j;
								cy[j] = p1;

								if(k == i - 1)
								{
									kflag = 0;
									inub = i + 1;
									s2_re = cy[j].re;
									s2_im = cy[j].im;
									j = -j;
									s1 = cy[j + 1];

									if(inub <= inu)
										goto_mw225 = true;
									else
									{
										s1.re = s2_re;
										s1.im = s2_im;
										goto_mw240 = true;
									}

									exitg1 = true;
								}
								else
								{
									k = i;
									b_guard1 = true;
								}
							}
							else
								b_guard2 = true;
						}
						else
							b_guard2 = true;

						if(b_guard2)
						{
							if(etest >= 350.46089684722295)
							{
								a2 -= 700.92179369444591;
								s1.re = cs_re * 3.9222272510438042E-305 - tm * 0.0;
								s1.im = cs_re * 0.0 + tm * 3.9222272510438042E-305;
								etest = s2_re;
								s2_re = s2_re * 3.9222272510438042E-305 - s2_im * 0.0;
								s2_im = etest * 0.0 + s2_im * 3.9222272510438042E-305;
								zd.re = a2;
								zd.im = yy;
							}

							b_guard1 = true;
						}

						if(b_guard1)
							i++;
					}

					if(goto_mw225 || goto_mw240)
						goto_mw110 = true;
					else
						goto_mw110 = false;

					if(!goto_mw110)
					{
						s1.re = s2_re;
						s1.im = s2_im;
						goto_mw270 = true;
					}
				}
				else
					goto_mw225 = true;
			}

			if(goto_mw225 || goto_mw240 || goto_mw270)
				goto_mw110 = true;
			else
				goto_mw110 = false;

			if(!goto_mw110)
			{
				s1.re = s2_re;
				s1.im = s2_im;
				zd = z;

				if(iflag == 1)
					;
				else
					goto_mw240 = true;
			}
		}
		else
			goto_mw225 = true;

		if(goto_mw225 || goto_mw240)
		{
			if(goto_mw225)
			{
				p1.re = dv7[kflag];
				rk = dv8[kflag];

				while(inub <= inu)
				{
					cs_re = s2_re;
					tm = s2_im;
					etest = s2_re;
					s2_re = (s2_re * ck_re - s2_im * ck_im) + s1.re;
					s2_im = (etest * ck_im + s2_im * ck_re) + s1.im;
					s1.re = cs_re;
					s1.im = tm;
					ck_re += rz_re;
					ck_im += rz_im;

					if(kflag + 1 < 3)
					{
						p2_re = s2_re * p1.re - s2_im * 0.0;
						p2_im = s2_re * 0.0 + s2_im * p1.re;
						etest = fabs(p2_re);
						a2 = fabs(p2_im);

						if((etest >= a2) || rtIsNaN(a2))
							b_etest = etest;
						else
							b_etest = a2;

						if(b_etest > rk)
						{
							kflag++;
							rk = dv8[kflag];
							s1.re = dv9[kflag] * (cs_re * p1.re - tm * 0.0);
							s1.im = dv9[kflag] * (cs_re * 0.0 + tm * p1.re);
							s2_re = dv9[kflag] * p2_re;
							s2_im = dv9[kflag] * p2_im;
							p1.re = dv7[kflag];
						}
					}

					inub++;
				}

				s1.re = s2_re;
				s1.im = s2_im;
			}

			y->re = dv7[kflag] * s1.re;
			y->im = dv7[kflag] * s1.im;
		}
		else
		{
			*nz = ckscl(zd, &s1, 700.92179369444591);
			*y = s1;

			if(1 <= *nz)
				;
			else
			{
				y->re = 2.2204460492503131E-16 * s1.re;
				y->im = 2.2204460492503131E-16 * s1.im;
				i1 = 1LL - *nz;

				if(i1 > 2147483647LL)
					i1 = 2147483647LL;
				else
				{
					if(i1 < -2147483648LL)
						i1 = -2147483648LL;
				}

				if((int)i1 >= 2)
				{
					y->re *= 2.2204460492503131E-16;
					y->im *= 2.2204460492503131E-16;
				}
			}
		}
	}
}

int cmlri(complex z, double fnu, int kode, int nin, complex *y, double tol)
{
	int nz;
	double az;
	double flooraz;
	int iaz;
	double fixfnu;
	int ifnu;
	int icounter;
	int inu;
	double ck_re;
	double tst;
	double ck_im;
	double ap;
	double ack;
	double rz_re;
	double rho;
	double rz_im;
	double p1_re;
	double p1_im;
	double p2_re;
	double p2_im;
	bool earlyExit;
	int i;
	bool exitg2;
	int kcounter;
	double pt_re;
	bool guard1 = false;
	double pt_im;
	int itime;
	bool exitg1;
	bool b_guard1 = false;
	nz = 0;
	az = rt_hypotd_snf(z.re, z.im);
	flooraz = floor(az);
	iaz = (int)flooraz;

	if(fnu < 0.0)
		fixfnu = ceil(fnu);
	else
		fixfnu = floor(fnu);

	ifnu = (int)fixfnu;

	if(nin <= 1)
		icounter = nin;
	else
		icounter = 1;

	inu = (ifnu + icounter) - 1;

	if(z.im == 0.0)
	{
		ck_re = (flooraz + 1.0) / z.re;
		ck_im = 0.0;
	}
	else if(z.re == 0.0)
	{
		ck_re = 0.0;
		ck_im = -((flooraz + 1.0) / z.im);
	}
	else
	{
		tst = fabs(z.re);
		ap = fabs(z.im);

		if(tst > ap)
		{
			ack = z.im / z.re;
			rho = z.re + ack * z.im;
			ck_re = ((flooraz + 1.0) + ack * 0.0) / rho;
			ck_im = (0.0 - ack * (flooraz + 1.0)) / rho;
		}
		else if(ap == tst)
		{
			if(z.re > 0.0)
				ack = 0.5;
			else
				ack = -0.5;

			if(z.im > 0.0)
				rho = 0.5;
			else
				rho = -0.5;

			ck_re = ((flooraz + 1.0) * ack + 0.0 * rho) / tst;
			ck_im = (0.0 * ack - (flooraz + 1.0) * rho) / tst;
		}
		else
		{
			ack = z.re / z.im;
			rho = z.im + ack * z.re;
			ck_re = ack * (flooraz + 1.0) / rho;
			ck_im = (ack * 0.0 - (flooraz + 1.0)) / rho;
		}
	}

	if(z.im == 0.0)
	{
		rz_re = 2.0 / z.re;
		rz_im = 0.0;
	}
	else if(z.re == 0.0)
	{
		rz_re = 0.0;
		rz_im = -(2.0 / z.im);
	}
	else
	{
		tst = fabs(z.re);
		ap = fabs(z.im);

		if(tst > ap)
		{
			ack = z.im / z.re;
			rho = z.re + ack * z.im;
			rz_re = (2.0 + ack * 0.0) / rho;
			rz_im = (0.0 - ack * 2.0) / rho;
		}
		else if(ap == tst)
		{
			if(z.re > 0.0)
				ack = 0.5;
			else
				ack = -0.5;

			if(z.im > 0.0)
				rho = 0.5;
			else
				rho = -0.5;

			rz_re = (2.0 * ack + 0.0 * rho) / tst;
			rz_im = (0.0 * ack - 2.0 * rho) / tst;
		}
		else
		{
			ack = z.re / z.im;
			rho = z.im + ack * z.re;
			rz_re = ack * 2.0 / rho;
			rz_im = (ack * 0.0 - 2.0) / rho;
		}
	}

	p1_re = 0.0;
	p1_im = 0.0;
	p2_re = 1.0;
	p2_im = 0.0;
	ack = ((flooraz + 1.0) + 1.0) / az;
	rho = ack + sqrt(ack * ack - 1.0);
	ack = rho * rho;
	tst = (ack + ack) / ((ack - 1.0) * (rho - 1.0)) / tol;
	flooraz++;
	earlyExit = true;
	icounter = 1;
	i = 1;
	exitg2 = false;

	while((!exitg2) && (i < 81))
	{
		icounter++;
		pt_re = p2_re;
		pt_im = p2_im;
		ack = ck_re * p2_im + ck_im * p2_re;
		p2_re = p1_re - (ck_re * p2_re - ck_im * p2_im);
		p2_im = p1_im - ack;
		p1_re = pt_re;
		p1_im = pt_im;
		ck_re += rz_re;
		ck_im += rz_im;

		if(rt_hypotd_snf(p2_re, p2_im) > tst * flooraz * flooraz)
		{
			earlyExit = false;
			exitg2 = true;
		}
		else
		{
			flooraz++;
			i++;
		}
	}

	if(earlyExit)
		nz = -2;
	else
	{
		kcounter = 1;
		guard1 = false;

		if(inu >= iaz)
		{
			p1_re = 0.0;
			p1_im = 0.0;
			p2_re = 1.0;
			p2_im = 0.0;

			if(z.im == 0.0)
			{
				ck_re = ((double)inu + 1.0) / z.re;
				ck_im = 0.0;
			}
			else if(z.re == 0.0)
			{
				if((double)inu + 1.0 == 0.0)
				{
					ck_re = 0.0 / z.im;
					ck_im = 0.0;
				}
				else
				{
					ck_re = 0.0;
					ck_im = -(((double)inu + 1.0) / z.im);
				}
			}
			else
			{
				tst = fabs(z.re);
				ap = fabs(z.im);

				if(tst > ap)
				{
					ack = z.im / z.re;
					rho = z.re + ack * z.im;
					ck_re = (((double)inu + 1.0) + ack * 0.0) / rho;
					ck_im = (0.0 - ack * ((double)inu + 1.0)) / rho;
				}
				else if(ap == tst)
				{
					if(z.re > 0.0)
						ack = 0.5;
					else
						ack = -0.5;

					if(z.im > 0.0)
						rho = 0.5;
					else
						rho = -0.5;

					ck_re = (((double)inu + 1.0) * ack + 0.0 * rho) / tst;
					ck_im = (0.0 * ack - ((double)inu + 1.0) * rho) / tst;
				}
				else
				{
					ack = z.re / z.im;
					rho = z.im + ack * z.re;
					ck_re = ack * ((double)inu + 1.0) / rho;
					ck_im = (ack * 0.0 - ((double)inu + 1.0)) / rho;
				}
			}

			tst = sqrt(((double)inu + 1.0) / az / tol);
			itime = 1;
			earlyExit = true;
			i = 1;
			exitg1 = false;

			while((!exitg1) && (i < 81))
			{
				kcounter++;
				pt_re = p2_re;
				pt_im = p2_im;
				ack = ck_re * p2_im + ck_im * p2_re;
				p2_re = p1_re - (ck_re * p2_re - ck_im * p2_im);
				p2_im = p1_im - ack;
				p1_re = pt_re;
				p1_im = pt_im;
				ck_re += rz_re;
				ck_im += rz_im;
				ap = rt_hypotd_snf(p2_re, p2_im);
				b_guard1 = false;

				if(ap >= tst * flooraz * flooraz)
				{
					if(itime == 2)
					{
						earlyExit = false;
						exitg1 = true;
					}
					else
					{
						ack = rt_hypotd_snf(ck_re, ck_im);
						ack += sqrt(ack * ack - 1.0);
						rho = ap / rt_hypotd_snf(pt_re, pt_im);

						if((ack <= rho) || rtIsNaN(rho))
							rho = ack;

						tst *= sqrt(rho / (rho * rho - 1.0));
						itime = 2;
						b_guard1 = true;
					}
				}
				else
					b_guard1 = true;

				if(b_guard1)
					i++;
			}

			if(earlyExit)
				nz = -2;
			else
				guard1 = true;
		}
		else
			guard1 = true;

		if(guard1)
		{
			itime = icounter + iaz;
			icounter = kcounter + inu;

			if(itime >= icounter)
				icounter = itime;

			az = icounter;
			p1_re = 0.0;
			p1_im = 0.0;
			p2_re = 2.2250738585072014E-305 / tol;
			p2_im = 0.0;
			fixfnu = fnu - fixfnu;
			flooraz = fixfnu + fixfnu;
			ack = ((double)icounter + flooraz) + 1.0;
			gammaln(&ack);
			rho = (double)icounter + 1.0;
			gammaln(&rho);
			ap = flooraz + 1.0;
			gammaln(&ap);
			ap = exp((ack - rho) - ap);
			ck_re = 0.0;
			ck_im = 0.0;
			icounter -= inu;

			for(i = 1; i <= icounter; i++)
			{
				pt_re = p2_re;
				pt_im = p2_im;
				ack = (az + fixfnu) * rz_re;
				rho = (az + fixfnu) * rz_im;
				tst = ack * p2_im + rho * p2_re;
				p2_re = p1_re + (ack * p2_re - rho * p2_im);
				p2_im = p1_im + tst;
				p1_re = pt_re;
				p1_im = pt_im;
				ack = ap * (1.0 - flooraz / (az + flooraz));
				ck_re += (ack + ap) * pt_re;
				ck_im += (ack + ap) * pt_im;
				ap = ack;
				az--;
			}

			y->re = p2_re;
			y->im = p2_im;

#if 0
			printf("y->re = %lf\n", y->re);
			printf("y->im = %lf\n", y->im);
#endif

			if(ifnu > 0)
			{
				for(i = 1; i <= ifnu; i++)
				{
					pt_re = p2_re;
					pt_im = p2_im;
					ack = (az + fixfnu) * rz_re;
					rho = (az + fixfnu) * rz_im;
					tst = ack * p2_im + rho * p2_re;
					p2_re = p1_re + (ack * p2_re - rho * p2_im);
					p2_im = p1_im + tst;
					p1_re = pt_re;
					p1_im = pt_im;
					ack = ap * (1.0 - flooraz / (az + flooraz));
					ck_re += (ack + ap) * pt_re;
					ck_im += (ack + ap) * pt_im;
					ap = ack;
					az--;
				}
			}

			pt_re = z.re;
			pt_im = z.im;

			if(kode == 2)
			{
				pt_re = z.re - z.re;
				pt_im = z.im;
			}

			if((rz_im == 0.0) && rtIsNaN(rz_re))
				printf("IsNaN!\n");
			else if((fabs(rz_re) > 8.9884656743115785E+307) || (fabs(rz_im) > 8.9884656743115785E+307))
			{
				printf("fabs(rz_re) > 8.9884656743115785E+307\n");
				ack = rz_re;
				rz_re = log(rt_hypotd_snf(rz_re / 2.0, rz_im / 2.0)) + 0.69314718055994529;
				rz_im = rt_atan2d_snf(rz_im, ack);
			}
			else
			{
				printf("else fabs(rz_re) > 8.9884656743115785E+307\n");
				ack = rz_re;
				rz_re = log(rt_hypotd_snf(rz_re, rz_im));
				rz_im = rt_atan2d_snf(rz_im, ack);
			}

			ap = 1.0 + fixfnu;
			gammaln(&ap);
			pt_re = ((-fixfnu * rz_re - -0.0 * rz_im) + pt_re) - ap;
			pt_im += -fixfnu * rz_im + -0.0 * rz_re;
			p2_re += ck_re;
			p2_im += ck_im;
			p1_re = 1.0 / rt_hypotd_snf(p2_re, p2_im);

			if(rtIsInf(pt_im) && rtIsInf(pt_re) && (pt_re < 0.0))
			{
				printf("IsInf!\n");
				rz_re = 0.0;
				rz_im = 0.0;
			}
			else
			{
				printf("else IsInf!\n");
				ack = exp(pt_re / 2.0);
				rz_re = ack * (ack * cos(pt_im));
				rz_im = ack * (ack * sin(pt_im));
			}

			ack = rz_re * p1_re - rz_im * 0.0;
			rz_im = rz_re * 0.0 + rz_im * p1_re;
			rho = p2_re * p1_re + p2_im * 0.0;
			p2_im = p2_re * 0.0 - p2_im * p1_re;
			rz_re = ack * rho - rz_im * p2_im;
			rz_im = ack * p2_im + rz_im * rho;
			ack = y->re;
			rho = y->im;
			y->re = ack * rz_re - rho * rz_im;
			y->im = ack * rz_im + rho * rz_re;

#if 1
			printf("y->re = %lf\n", y->re);
			printf("y->im = %lf\n", y->im);
#endif
		}
	}

	return nz;
}

void b_cosh(complex *x)
{
	double x_re;
	double x_im;

	if(x->im == 0.0)
	{
		x->re = cosh(x->re);
		x->im = 0.0;
	}
	else
	{
		x_re = x->re;
		x_im = x->im;
		x->re = cosh(x->re) * cos(x->im);
		x->im = sinh(x_re) * sin(x_im);
	}
}

void b_sinh(complex *x)
{
	double x_re;
	double x_im;

	if(x->im == 0.0)
	{
		x->re = sinh(x->re);
		x->im = 0.0;
	}
	else
	{
		x_re = x->re;
		x_im = x->im;
		x->re = sinh(x->re) * cos(x->im);
		x->im = cosh(x_re) * sin(x_im);
	}
}

int cseri(complex z, double fnu, int kode, int nin, complex *y, double tol, double elim, double alim)
{
	int nz;
	int n;
	double az;
	double crsc_re;
	bool iflag;
	double hz_re;
	double hz_im;
	double cz_re;
	double cz_im;
	double acz;
	complex ck;
	double ak;
	double ascle;
	double aa;
	double coef_re;
	double coef_im;
	double b_atol;
	int i;
	bool exitg1;
	double s1_re;
	double s1_im;

	if(nin <= 1)
		n = nin;
	else
		n = 1;

	nz = 0;
	az = rt_hypotd_snf(z.re, z.im);

	if(az == 0.0)
	{
		if(fnu == 0.0)
		{
			y->re = 1.0;
			y->im = 0.0;
		}
		else
		{
			y->re = 0.0;
			y->im = 0.0;
		}
	}
	else
	{
		crsc_re = 1.0;
		iflag = false;

		if (az < 2.2250738585072014E-305)
		{
			nz = n;

			if(fnu == 0.0)
				nz = n - 1;

			if(fnu == 0.0)
			{
				y->re = 1.0;
				y->im = 0.0;
			}
			else
			{
				y->re = 0.0;
				y->im = 0.0;
			}
		}
		else
		{
			hz_re = 0.5 * z.re;
			hz_im = 0.5 * z.im;

			if(az > 4.7170688552396617E-153)
			{
				cz_re = hz_re * hz_re - hz_im * hz_im;
				cz_im = hz_re * hz_im + hz_im * hz_re;
				acz = rt_hypotd_snf(cz_re, cz_im);
			}
			else
			{
				cz_re = 0.0;
				cz_im = 0.0;
				acz = 0.0;
			}

			ck.re = hz_re;
			ck.im = hz_im;

			if((hz_im == 0.0) && rtIsNaN(hz_re))
				printf("IsNaN detection\n");
			else if((fabs(hz_re) > 8.9884656743115785E+307) || (fabs(hz_im) > 8.9884656743115785E+307))
			{
				ck.re = log(rt_hypotd_snf(hz_re / 2.0, hz_im / 2.0)) + 0.69314718055994529;
				ck.im = rt_atan2d_snf(hz_im, hz_re);
			}
			else
			{
				ck.re = log(rt_hypotd_snf(hz_re, hz_im));
				ck.im = rt_atan2d_snf(hz_im, hz_re);
			}

			az = (fnu + (double)n) - 1.0;
			ak = az + 1.0;
			gammaln(&ak);
			ck.re = ck.re * az - ak;
			ck.im *= az;

			if(kode == 2)
				ck.re -= z.re;

			if(ck.re > -elim)
			{
				printf("elim\n");
				ascle = 0.0;

				if(ck.re <= -alim)
				{
					printf("alim\n");
					iflag = true;
					crsc_re = tol;
					ascle = 2.2250738585072014E-305 / tol;
				}

				aa = exp(ck.re);

				if(iflag)
					aa /= tol;

				coef_re = aa * cos(ck.im);
				coef_im = aa * sin(ck.im);
				b_atol = tol * acz / (az + 1.0);
				i = 1;
				exitg1 = false;

				while((!exitg1) && (i <= n))
				{
					az = fnu;
					s1_re = 1.0;
					s1_im = 0.0;

					if(!(acz < tol * (az + 1.0)))
					{
						ck.re = 1.0;
						ck.im = 0.0;
						ak = (az + 1.0) + 2.0;
						az++;
						aa = 2.0;

						do
						{
							hz_re = 1.0 / az;
							hz_im = ck.re;
							ck.re = hz_re * (ck.re * cz_re - ck.im * cz_im);
							ck.im = hz_re * (hz_im * cz_im + ck.im * cz_re);
							s1_re += ck.re;
							s1_im += ck.im;
							az += ak;
							ak += 2.0;
							aa = aa * acz * hz_re;
						}
						while(!!(aa > b_atol));
					}

					ck.re = s1_re * coef_re - s1_im * coef_im;
					ck.im = s1_re * coef_im + s1_im * coef_re;

					if(iflag && (cuchk(ck, ascle, tol) != 0))
						exitg1 = true;
					else
					{
						y->re = ck.re * crsc_re - ck.im * 0.0;
						y->im = ck.re * 0.0 + ck.im * crsc_re;
						i = 2;
					}
				}
			}
			else
			{
				nz = 1;
				y->re = 0.0;
				y->im = 0.0;

				if(acz > az)
					nz = -1;
			}
		}
	}

	return nz;
}

int cuchk(complex y, double ascle, double tol)
{
	int nz;
	double yr;
	double yi;
	double smallpart;
	yr = fabs(y.re);
	yi = fabs(y.im);

	if(yr > yi)
	{
		smallpart = yi;
		yi = yr;
	}
	else
		smallpart = yr;

	if((smallpart <= ascle) && (yi < smallpart / tol))
		nz = 1;
	else
		nz = 0;

	return nz;
}

double rt_powd_snf(double u0, double u1)
{
	double y;
	double d0;
	double d1;

	if(rtIsNaN(u0) || rtIsNaN(u1))
		y = rtNaN;
	else
	{
		d0 = fabs(u0);
		d1 = fabs(u1);

		if(rtIsInf(u1))
		{
			if(d0 == 1.0)
				y = rtNaN;
			else if (d0 > 1.0)
			{
				if(u1 > 0.0)
					y = rtInf;
				else
					y = 0.0;
			}
			else if(u1 > 0.0)
				y = 0.0;
			else
				y = rtInf;
		}
		else if(d1 == 0.0)
			y = 1.0;
		else if(d1 == 1.0)
		{
			if(u1 > 0.0)
				y = u0;
			else
				y = 1.0 / u0;
		}
		else if(u1 == 2.0)
			y = u0 * u0;
		else if((u1 == 0.5) && (u0 >= 0.0))
			y = sqrt(u0);
		else if((u0 < 0.0) && (u1 > floor(u1)))
			y = rtNaN;
		else
			y = pow(u0, u1);
	}

	return y;
}

void b_cunhj(complex z, double fnu, int ipmtr, double tol, complex *phi, complex *arg, complex *zeta1, complex *zeta2, complex *asum, complex *bsum)
{
	double rfnu;
	complex up[14];
	double ac;
	complex zb;
	double rfnu2;
	double fn23;
	double rfn13_re;
	complex w2;
	int k;
	complex w;
	complex p[30];
	double ap[30];
	double pp;
	complex suma;
	int l1;
	double brm;
	complex zc;
	bool exitg4;
	double btol;
	double d;
	static const double GAMA[30] = { 0.6299605249474366, 0.25198420997897464,
		0.15479030041565583, 0.11071306241615901, 0.085730939552739485,
		0.069716131695868433, 0.058608567189371359, 0.050469887353631067,
		0.044260058068915482, 0.039372066154350994, 0.035428319592445537,
		0.032181885750209825, 0.029464624079115768, 0.027158167711293448,
		0.025176827297386177, 0.023457075530607888, 0.021950839013490719,
		0.020621082823564625, 0.019438824089788084, 0.018381063380068317,
		0.017429321323196318, 0.016568583778661234, 0.015786528598791844,
		0.015072950149409559, 0.014419325083995464, 0.013818480573534178,
		0.013264337899427657, 0.012751712197049864, 0.012276154531876277,
		0.011833826239848241 };

	double zth_re;
	double zth_im;
	static const double BETA[210] = { 0.017998872141355329, 0.0055996491106438812,
		0.0028850140223113277, 0.0018009660676105393, 0.0012475311058919921,
		0.00092287887657293828, 0.00071443042172728737, 0.00057178728178970488,
		0.00046943100760648155, 0.00039323283546291665, 0.00033481888931829768,
		0.00028895214849575154, 0.0002522116155495733, 0.00022228058079888332,
		0.00019754183803306251, 0.00017683685501971802, 0.00015931689966182109,
		0.00014434793019733397, 0.00013144806811996539, 0.00012024544494930288,
		0.00011044914450459939, 0.00010182877074056726, 9.4199822420423752E-5,
		8.7413054575383449E-5, 8.1346626216280142E-5, 7.590022696462193E-5,
		7.0990630063415351E-5, 6.6548287484246817E-5, 6.25146958969275E-5,
		5.8840339442625178E-5, -0.0014928295321342917, -0.00087820470954638936,
		-0.00050291654957203464, -0.000294822138512746, -0.00017546399697078284,
		-0.00010400855046081644, -5.961419530464579E-5, -3.1203892907609836E-5,
		-1.2608973598023005E-5, -2.4289260857573037E-7, 8.059961654142736E-6,
		1.3650700926214739E-5, 1.7396412547292627E-5, 1.9867297884213378E-5,
		2.1446326379082263E-5, 2.2395465923245652E-5, 2.2896778381471263E-5,
		2.3078538981117782E-5, 2.3032197608090914E-5, 2.2823607372034874E-5,
		2.2500588110529241E-5, 2.2098101536199144E-5, 2.164184274481039E-5,
		2.1150764925622083E-5, 2.0638874978217072E-5, 2.0116524199708165E-5,
		1.9591345014117925E-5, 1.9068936791043675E-5, 1.8553371964163667E-5,
		1.8047572225967421E-5, 0.0005522130767212928, 0.00044793258155238465,
		0.00027952065399202059, 0.0001524681561984466, 6.932711056570436E-5,
		1.7625868306999139E-5, -1.3574499634326914E-5, -3.1797241335042717E-5,
		-4.188618616966934E-5, -4.6900488937914104E-5, -4.8766544741378735E-5,
		-4.8701003118673505E-5, -4.7475562089008661E-5, -4.5581305813862843E-5,
		-4.33309644511266E-5, -4.0923019315775034E-5, -3.848226386032213E-5,
		-3.6085716753541052E-5, -3.3779330612336739E-5, -3.1588856077210961E-5,
		-2.9526956175080731E-5, -2.7597891482833575E-5, -2.5800617466688372E-5,
		-2.4130835676128019E-5, -2.2582350951834605E-5, -2.1147965676891298E-5,
		-1.9820063888529493E-5, -1.8590987080106508E-5, -1.7453269984421023E-5,
		-1.63997823854498E-5, -0.0004746177965599598, -0.0004778645671473215,
		-0.00032039022806703763, -0.00016110501611996228, -4.2577810128543523E-5,
		3.4457129429496748E-5, 7.97092684075675E-5, 0.00010313823670827221,
		0.00011246677526220416, 0.00011310364210848139, 0.00010865163484877427,
		0.00010143795159766197, 9.29298396593364E-5, 8.4029313301609E-5,
		7.52727991349134E-5, 6.696325219757309E-5, 5.925645473231947E-5,
		5.2216930882697554E-5, 4.5853948516536063E-5, 4.0144551389148682E-5,
		3.5048173003132809E-5, 3.0515799503434667E-5, 2.6495611995051603E-5,
		2.2936363369099816E-5, 1.9789305666402162E-5, 1.7009198463641262E-5,
		1.45547428261524E-5, 1.2388664099587841E-5, 1.0477587607658323E-5,
		8.7917995497847932E-6, 0.00073646581057257841, 0.000872790805146194,
		0.00062261486257313506, 0.00028599815419430417, 3.8473767287936606E-6,
		-0.00018790600363697156, -0.00029760364659455455, -0.00034599812683265633,
		-0.00035338247091603773, -0.00033571563577504876, -0.00030432112478903981,
		-0.00026672272304761283, -0.00022765421412281953, -0.00018992261185456235,
		-0.00015505891859909386, -0.00012377824076187363, -9.6292614771764412E-5,
		-7.2517832771442527E-5, -5.2207002889563382E-5, -3.5034775051190054E-5,
		-2.0648976103555174E-5, -8.7010609684976711E-6, 1.136986866751003E-6,
		9.1642647412277879E-6, 1.5647778542887261E-5, 2.0822362948246685E-5,
		2.4892338100459516E-5, 2.8034050957414632E-5, 3.0398777462986191E-5,
		3.2115673140670063E-5, -0.0018018219196388571, -0.0024340296293804253,
		-0.001834226635498568, -0.00076220459635400974, 0.00023907947525692722,
		0.00094926611717688109, 0.0013446744970154036, 0.0014845749525944918,
		0.0014473233983061759, 0.0013026826128565718, 0.0011035159737564268,
		0.00088604744041979172, 0.00067307320816566542, 0.00047760387285658237,
		0.00030599192635878935, 0.00016031569459472162, 4.0074955527061327E-5,
		-5.6660746163525162E-5, -0.00013250618677298264, -0.00019029618798961406,
		-0.00023281145037693741, -0.00026262881146466884, -0.00028205046986759866,
		-0.00029308156319286116, -0.00029743596217631662, -0.00029655733423934809,
		-0.00029164736331209088, -0.00028369620383773418, -0.00027351231709567335,
		-0.00026175015580676858, 0.0063858589121205088, 0.00962374215806378,
		0.0076187806120700105, 0.0028321905554562804, -0.0020984135201272008,
		-0.0057382676421662646, -0.0077080424449541465, -0.0082101169226484437,
		-0.0076582452034690543, -0.006472097293910452, -0.0049913241200496648,
		-0.0034561228971313326, -0.0020178558001417079, -0.00075943068678196145,
		0.00028417363152385912, 0.0011089166758633741, 0.0017290149387272878,
		0.0021681259080268472, 0.0024535771049453972, 0.0026128182105833488,
		0.0026714103965627691, 0.0026520307339598045, 0.0025741165287728731,
		0.0024538912623609443, 0.002304600580717955, 0.0021368483768671267,
		0.0019589652847887091, 0.0017773700867945441, 0.0015969028076583906,
		0.0014211197566443854 };

	int l2;
	int ias;
	int ibs;
	double rtzta_re;
	int lr;
	double rtzta_im;
	bool exitg1;
	int ks;
	bool exitg3;
	bool exitg2;
	static const double ALFA[180] = { -0.0044444444444444444,
		-0.000922077922077922, -8.8489288489288488E-5, 0.00016592768783244973,
		0.00024669137274179289, 0.00026599558934625478, 0.00026182429706150096,
		0.00024873043734465562, 0.00023272104008323209, 0.00021636248571236508,
		0.00020073885876275234, 0.00018626763663754517, 0.0001730607759178765,
		0.00016109170592901574, 0.00015027477416090814, 0.00014050349739126979,
		0.0001316688165459228, 0.00012366744559825325, 0.00011640527147473791,
		0.00010979829837271337, 0.00010377241042299283, 9.8262607836936344E-5,
		9.321205172495032E-5, 8.857108524787117E-5, 8.4296310571570029E-5,
		8.0349754840779115E-5, 7.6698134535920737E-5, 7.3312215748177779E-5,
		7.0166262516314139E-5, 6.7237563379016026E-5, 0.000693735541354589,
		0.00023224174518292166, -1.419862735566912E-5, -0.00011644493167204864,
		-0.00015080355805304876, -0.00015512192491809622, -0.00014680975664646556,
		-0.00013381550386749137, -0.00011974497568425405, -0.00010618431920797402,
		-9.3769954989119444E-5, -8.2692304558819327E-5, -7.2937434815522126E-5,
		-6.4404235772101633E-5, -5.69611566009369E-5, -5.0473104430356164E-5,
		-4.4813486800888282E-5, -3.9868872771759884E-5, -3.5540053297204251E-5,
		-3.174142566090225E-5, -2.839967939041748E-5, -2.5452272063487058E-5,
		-2.2845929716472455E-5, -2.0535275310648061E-5, -1.848162176276661E-5,
		-1.6651933002139381E-5, -1.5017941298011949E-5, -1.3555403137904052E-5,
		-1.2243474647385812E-5, -1.1064188481130817E-5, -0.00035421197145774384,
		-0.00015616126394515941, 3.0446550359493642E-5, 0.00013019865577324269,
		0.00016747110669971228, 0.00017022258768359256, 0.00015650142760859472,
		0.00013633917097744512, 0.00011488669202982512, 9.4586909303468817E-5,
		7.6449841925089825E-5, 6.0757033496519734E-5, 4.7439429929050881E-5,
		3.6275751200534429E-5, 2.6993971497922491E-5, 1.9321093824793926E-5,
		1.3005667479396321E-5, 7.8262086674449658E-6, 3.5925748581935159E-6,
		1.4404004981425182E-7, -2.6539676969793912E-6, -4.9134686709848593E-6,
		-6.7273929609124832E-6, -8.17269379678658E-6, -9.313047150935612E-6,
		-1.0201141879801643E-5, -1.0880596251059288E-5, -1.1387548150960355E-5,
		-1.1751967567455642E-5, -1.1998736487094414E-5, 0.00037819419920177291,
		0.00020247195276181616, -6.3793850631886236E-5, -0.0002385982306030059,
		-0.00031091625602736159, -0.00031368011524757634, -0.00027895027379132341,
		-0.00022856408261914138, -0.00017524528034084676, -0.00012554406306069035,
		-8.2298287282020835E-5, -4.6286073058811649E-5, -1.7233430236696227E-5,
		5.6069048230460226E-6, 2.313954431482868E-5, 3.6264274585679393E-5,
		4.5800612449018877E-5, 5.2459529495911405E-5, 5.6839620854581527E-5,
		5.9434982039310406E-5, 6.0647852757842175E-5, 6.0802390778843649E-5,
		6.0157789453946036E-5, 5.8919965734469847E-5, 5.72515823777593E-5,
		5.5280437558585257E-5, 5.3106377380288019E-5, 5.080693020123257E-5,
		4.8441864762009484E-5, 4.6056858160747536E-5, -0.00069114139728829421,
		-0.00042997663305887192, 0.000183067735980039, 0.00066008814754201417,
		0.0008759649699511859, 0.00087733523595823551, 0.00074936958537899067,
		0.000563832329756981, 0.00036805931997144317, 0.00018846453551445559,
		3.7066305766490415E-5, -8.28520220232137E-5, -0.000172751952869173,
		-0.00023631487360587297, -0.00027796615069490668, -0.00030207951415545694,
		-0.00031259471264382012, -0.00031287255875806717, -0.0003056780384663244,
		-0.00029322647061455731, -0.0002772556555829348, -0.00025910392846703172,
		-0.00023978401439648034, -0.00022004826004542284, -0.00020044391109497149,
		-0.00018135869221097068, -0.00016305767447865748, -0.00014571267217520584,
		-0.0001294254219839246, -0.00011424569194244596, 0.0019282196424877589,
		0.0013559257630202223, -0.000717858090421303, -0.0025808480257527035,
		-0.0034927113082616847, -0.0034698629934096061, -0.0028228523335131019,
		-0.0018810307640489134, -0.00088953171838394764, 3.8791210263103525E-6,
		0.00072868854011969139, 0.0012656637305345775, 0.0016251815837267443,
		0.0018320315321637317, 0.0019158838899052792, 0.0019058884675554615,
		0.0018279898242182574, 0.0017038950642112153, 0.0015509712717109768,
		0.0013826142185227616, 0.0012088142423006478, 0.0010367653263834496,
		0.00087143791806861914, 0.000716080155297701, 0.00057263700255812935,
		0.00044208981946580229, 0.00032472494850309055, 0.00022034204273024659,
		0.00012841289840135388, 4.8200592455209545E-5 };

	double tfn_im;
	double rzth_re;
	double rzth_im;
	int kp1;
	int l;
	complex cr[14];
	complex dr[14];
	bool exitg5;
	static const double C[105] = { 1.0, -0.20833333333333334, 0.125,
		0.3342013888888889, -0.40104166666666669, 0.0703125, -1.0258125964506173,
		1.8464626736111112, -0.8912109375, 0.0732421875, 4.6695844234262474,
		-11.207002616222994, 8.78912353515625, -2.3640869140625, 0.112152099609375,
		-28.212072558200244, 84.636217674600729, -91.818241543240021,
		42.534998745388457, -7.3687943594796321, 0.22710800170898438,
		212.57013003921713, -765.25246814118168, 1059.9904525279999,
		-699.57962737613252, 218.19051174421159, -26.491430486951554,
		0.57250142097473145, -1919.4576623184071, 8061.7221817373093,
		-13586.550006434138, 11655.393336864534, -5305.646978613403,
		1200.9029132163525, -108.09091978839466, 1.7277275025844574,
		20204.291330966149, -96980.598388637518, 192547.00123253153,
		-203400.17728041555, 122200.46498301746, -41192.65496889755,
		7109.5143024893641, -493.915304773088, 6.074042001273483,
		-242919.18790055133, 1.3117636146629772E+6, -2.9980159185381066E+6,
		3.7632712976564039E+6, -2.8135632265865342E+6, 1.2683652733216248E+6,
		-331645.17248456361, 45218.768981362729, -2499.8304818112097,
		24.380529699556064, 3.2844698530720379E+6, -1.9706819118432228E+7,
		5.0952602492664643E+7, -7.4105148211532652E+7, 6.6344512274729028E+7,
		-3.7567176660763353E+7, 1.3288767166421818E+7, -2.7856181280864547E+6,
		308186.40461266239, -13886.08975371704, 110.01714026924674,
		-4.932925366450996E+7, 3.2557307418576574E+8, -9.394623596815784E+8,
		1.55359689957058E+9, -1.6210805521083372E+9, 1.1068428168230145E+9,
		-4.9588978427503031E+8, 1.4206290779753309E+8, -2.447406272573873E+7,
		2.2437681779224495E+6, -84005.433603024081, 551.33589612202059,
		8.1478909611831212E+8, -5.8664814920518475E+9, 1.8688207509295826E+10,
		-3.4632043388158775E+10, 4.1280185579753975E+10, -3.3026599749800724E+10,
		1.79542137311556E+10, -6.5632937926192846E+9, 1.5592798648792574E+9,
		-2.2510566188941526E+8, 1.7395107553978164E+7, -549842.32757228869,
		3038.0905109223841, -1.4679261247695616E+10, 1.144982377320258E+11,
		-3.9909617522446649E+11, 8.1921866954857727E+11, -1.0983751560812233E+12,
		1.0081581068653821E+12, -6.4536486924537646E+11, 2.8790064990615057E+11,
		-8.786707217802327E+10, 1.7634730606834969E+10, -2.1671649832237949E+9,
		1.4315787671888897E+8, -3.8718334425726128E+6, 18257.755474293175 };

	static const double BR[14] = { 1.0, -0.14583333333333334,
		-0.098741319444444448, -0.14331205391589505, -0.31722720267841353,
		-0.94242914795712029, -3.5112030408263544, -15.727263620368046,
		-82.281439097185938, -492.3553705236705, -3316.2185685479726,
		-24827.674245208589, -204526.58731512979, -1.83844491706821E+6 };

	static const double AR[14] = { 1.0, 0.10416666666666667, 0.083550347222222224,
		0.12822657455632716, 0.29184902646414046, 0.88162726744375763,
		3.3214082818627677, 14.995762986862555, 78.923013011586519,
		474.45153886826432, 3207.490090890662, 24086.549640874004, 198923.1191695098,
		1.7919020077753437E+6 };

	rfnu = 1.0 / fnu;
	memset(&up[0], 0, 14U * sizeof(complex));
	ac = fnu * 2.2250738585072014E-305;
	asum->re = 0.0;
	asum->im = 0.0;
	bsum->re = 0.0;
	bsum->im = 0.0;

	if((fabs(z.re) <= ac) && (fabs(z.im) <= ac))
	{
		zeta1->re = 1402.9773265065639 + fnu;
		zeta1->im = 0.0;
		zeta2->re = fnu;
		zeta2->im = 0.0;
		phi->re = 1.0;
		phi->im = 0.0;
		arg->re = 1.0;
		arg->im = 0.0;
	}
	else
	{
		zb.re = rfnu * z.re;
		zb.im = rfnu * z.im;
		rfnu2 = rfnu * rfnu;
		ac = rt_powd_snf(fnu, 0.33333333333333331);
		fn23 = ac * ac;
		rfn13_re = 1.0 / ac;
		w2.re = 1.0 - (zb.re * zb.re - zb.im * zb.im);
		w2.im = 0.0 - (zb.re * zb.im + zb.im * zb.re);
		ac = rt_hypotd_snf(w2.re, w2.im);

		if(ac > 0.25)
		{
			w = w2;
			b_sqrt(&w);
			ac = w.re;
			pp = w.im;

			if(w.re < 0.0)
				ac = 0.0;

			if (w.im < 0.0)
				pp = 0.0;

			w.re = ac;
			w.im = pp;

			if(zb.im == 0.0)
			{
				if(pp == 0.0)
				{
					zc.re = (1.0 + ac) / zb.re;
					zc.im = 0.0;
				}
				else if(1.0 + ac == 0.0)
				{
					zc.re = 0.0;
					zc.im = pp / zb.re;
				}
				else
				{
					zc.re = (1.0 + ac) / zb.re;
					zc.im = pp / zb.re;
				}
			}
			else if(zb.re == 0.0)
			{
				if(1.0 + ac == 0.0)
				{
					zc.re = pp / zb.im;
					zc.im = 0.0;
				}
				else if (pp == 0.0)
				{
					zc.re = 0.0;
					zc.im = -((1.0 + ac) / zb.im);
				}
				else
				{
					zc.re = pp / zb.im;
					zc.im = -((1.0 + ac) / zb.im);
				}
			}
			else
			{
				brm = fabs(zb.re);
				btol = fabs(zb.im);

				if(brm > btol)
				{
					btol = zb.im / zb.re;
					d = zb.re + btol * zb.im;
					zc.re = ((1.0 + ac) + btol * pp) / d;
					zc.im = (pp - btol * (1.0 + ac)) / d;
				}
				else if(btol == brm)
				{
					if(zb.re > 0.0)
						btol = 0.5;
					else
						btol = -0.5;

					if(zb.im > 0.0)
						d = 0.5;
					else
						d = -0.5;

					zc.re = ((1.0 + ac) * btol + pp * d) / brm;
					zc.im = (pp * btol - (1.0 + ac) * d) / brm;
				}
				else
				{
					btol = zb.re / zb.im;
					d = zb.im + btol * zb.re;
					zc.re = (btol * (1.0 + ac) + pp) / d;
					zc.im = (btol * pp - (1.0 + ac)) / d;
				}
			}

			b_log(&zc);
			ac = zc.re;
			pp = zc.im;

			if(zc.re < 0.0)
				ac = 0.0;

			if(zc.im < 0.0)
				pp = 0.0;
			else
			{
				if(zc.im > 1.5707963267948966)
					pp = 1.5707963267948966;
			}

			zth_re = 1.5 * (ac - w.re);
			zth_im = 1.5 * (pp - w.im);
			zeta1->re = fnu * ac;
			zeta1->im = fnu * pp;
			zeta2->re = fnu * w.re;
			zeta2->im = fnu * w.im;

			if((zth_re >= 0.0) && (zth_im < 0.0))
				ac = 4.71238898038469;
			else if (zth_re != 0.0)
			{
				ac = atan(zth_im / zth_re);

				if(zth_re < 0.0)
					ac += 3.1415926535897931;
			}
			else
				ac = 1.5707963267948966;

			pp = rt_powd_snf(rt_hypotd_snf(zth_re, zth_im), 0.66666666666666663);
			ac *= 0.66666666666666663;
			zb.re = pp * cos(ac);
			zb.im = pp * sin(ac);

			if(zb.im < 0.0)
				zb.im = 0.0;

			arg->re = fn23 * zb.re;
			arg->im = fn23 * zb.im;

			if(zb.im == 0.0)
			{
				if(zth_im == 0.0)
				{
					rtzta_re = zth_re / zb.re;
					rtzta_im = 0.0;
				}
				else if(zth_re == 0.0)
				{
					rtzta_re = 0.0;
					rtzta_im = zth_im / zb.re;
				}
				else
				{
					rtzta_re = zth_re / zb.re;
					rtzta_im = zth_im / zb.re;
				}
			}
			else if(zb.re == 0.0)
			{
				if(zth_re == 0.0)
				{
					rtzta_re = zth_im / zb.im;
					rtzta_im = 0.0;
				}
				else if(zth_im == 0.0)
				{
					rtzta_re = 0.0;
					rtzta_im = -(zth_re / zb.im);
				}
				else
				{
					rtzta_re = zth_im / zb.im;
					rtzta_im = -(zth_re / zb.im);
				}
			}
			else
			{
				brm = fabs(zb.re);
				btol = zb.im;

				if(brm > btol)
				{
					btol = zb.im / zb.re;
					d = zb.re + btol * zb.im;
					rtzta_re = (zth_re + btol * zth_im) / d;
					rtzta_im = (zth_im - btol * zth_re) / d;
				}
				else if (btol == brm)
				{
					if(zb.re > 0.0)
						btol = 0.5;
					else
						btol = -0.5;

					if(zb.im > 0.0)
						d = 0.5;
					else
						d = -0.5;

					rtzta_re = (zth_re * btol + zth_im * d) / brm;
					rtzta_im = (zth_im * btol - zth_re * d) / brm;
				}
				else
				{
					btol = zb.re / zb.im;
					d = zb.im + btol * zb.re;
					rtzta_re = (btol * zth_re + zth_im) / d;
					rtzta_im = (btol * zth_im - zth_re) / d;
				}
			}

			if(w.im == 0.0)
			{
				if(rtzta_im == 0.0)
				{
					suma.re = rtzta_re / w.re;
					suma.im = 0.0;
				}
				else if(rtzta_re == 0.0)
				{
					suma.re = 0.0;
					suma.im = rtzta_im / w.re;
				}
				else
				{
					suma.re = rtzta_re / w.re;
					suma.im = rtzta_im / w.re;
				}
			}
			else if(w.re == 0.0)
			{
				if(rtzta_re == 0.0)
				{
					suma.re = rtzta_im / w.im;
					suma.im = 0.0;
				}
				else if(rtzta_im == 0.0)
				{
					suma.re = 0.0;
					suma.im = -(rtzta_re / w.im);
				}
				else
				{
					suma.re = rtzta_im / w.im;
					suma.im = -(rtzta_re / w.im);
				}
			}
			else
			{
				brm = fabs(w.re);
				btol = fabs(w.im);

				if(brm > btol)
				{
					btol = w.im / w.re;
					d = w.re + btol * w.im;
					suma.re = (rtzta_re + btol * rtzta_im) / d;
					suma.im = (rtzta_im - btol * rtzta_re) / d;
				}
				else if(btol == brm)
				{
					if(w.re > 0.0)
						btol = 0.5;
					else
						btol = -0.5;

					if(w.im > 0.0)
						d = 0.5;
					else
						d = -0.5;

					suma.re = (rtzta_re * btol + rtzta_im * d) / brm;
					suma.im = (rtzta_im * btol - rtzta_re * d) / brm;
				}
				else
				{
					btol = w.re / w.im;
					d = w.im + btol * w.re;
					suma.re = (btol * rtzta_re + rtzta_im) / d;
					suma.im = (btol * rtzta_im - rtzta_re) / d;
				}
			}

			zb.re = suma.re + suma.re;
			zb.im = suma.im + suma.im;
			b_sqrt(&zb);
			phi->re = zb.re * rfn13_re - zb.im * 0.0;
			phi->im = zb.re * 0.0 + zb.im * rfn13_re;

			if (ipmtr == 1)
				;
			else
			{
				if(w.im == 0.0)
				{
					fn23 = rfnu / w.re;
					tfn_im = 0.0;
				}
				else if(w.re == 0.0)
				{
					if(rfnu == 0.0)
					{
						fn23 = 0.0 / w.im;
						tfn_im = 0.0;
					}
					else
					{
						fn23 = 0.0;
						tfn_im = -(rfnu / w.im);
					}
				}
				else
				{
					brm = fabs(w.re);
					btol = fabs(w.im);

					if(brm > btol)
					{
						btol = w.im / w.re;
						d = w.re + btol * w.im;
						fn23 = (rfnu + btol * 0.0) / d;
						tfn_im = (0.0 - btol * rfnu) / d;
					}
					else if(btol == brm)
					{
						if(w.re > 0.0)
							btol = 0.5;
						else
							btol = -0.5;

						if(w.im > 0.0)
							d = 0.5;
						else
							d = -0.5;

						fn23 = (rfnu * btol + 0.0 * d) / brm;
						tfn_im = (0.0 * btol - rfnu * d) / brm;
					}
					else
					{
						btol = w.re / w.im;
						d = w.im + btol * w.re;
						fn23 = btol * rfnu / d;
						tfn_im = (btol * 0.0 - rfnu) / d;
					}
				}

				if(zth_im == 0.0)
				{
					rzth_re = rfnu / zth_re;
					rzth_im = 0.0;
				}
				else if(zth_re == 0.0)
				{
					if(rfnu == 0.0)
					{
						rzth_re = 0.0 / zth_im;
						rzth_im = 0.0;
					}
					else
					{
						rzth_re = 0.0;
						rzth_im = -(rfnu / zth_im);
					}
				}
				else
				{
					brm = fabs(zth_re);
					btol = fabs(zth_im);

					if(brm > btol)
					{
						btol = zth_im / zth_re;
						d = zth_re + btol * zth_im;
						rzth_re = (rfnu + btol * 0.0) / d;
						rzth_im = (0.0 - btol * rfnu) / d;
					}
					else if(btol == brm)
					{
						if(zth_re > 0.0)
							btol = 0.5;
						else
							btol = -0.5;

						if(zth_im > 0.0)
							d = 0.5;
						else
							d = -0.5;

						rzth_re = (rfnu * btol + 0.0 * d) / brm;
						rzth_im = (0.0 * btol - rfnu * d) / brm;
					}
					else
					{
						btol = zth_re / zth_im;
						d = zth_im + btol * zth_re;
						rzth_re = btol * rfnu / d;
						rzth_im = (btol * 0.0 - rfnu) / d;
					}
				}

				zc.re = 0.10416666666666667 * rzth_re;
				zc.im = 0.10416666666666667 * rzth_im;

				if(w2.im == 0.0)
				{
					w.re = 1.0 / w2.re;
					w.im = 0.0;
				}
				else if(w2.re == 0.0)
				{
					w.re = 0.0;
					w.im = -(1.0 / w2.im);
				}
				else
				{
					brm = fabs(w2.re);
					btol = fabs(w2.im);

					if(brm > btol)
					{
						btol = w2.im / w2.re;
						d = w2.re + btol * w2.im;
						w.re = (1.0 + btol * 0.0) / d;
						w.im = (0.0 - btol) / d;
					}
					else if(btol == brm)
					{
						if(w2.re > 0.0)
							btol = 0.5;
						else
							btol = -0.5;

						if(w2.im > 0.0)
							d = 0.5;
						else
							d = -0.5;

						w.re = (btol + 0.0 * d) / brm;
						w.im = (0.0 * btol - d) / brm;
					}
					else
					{
						btol = w2.re / w2.im;
						d = w2.im + btol * w2.re;
						w.re = btol / d;
						w.im = (btol * 0.0 - 1.0) / d;
					}
				}

				ac = w.re * -0.20833333333333334 + 0.125;
				pp = w.im * -0.20833333333333334;
				up[1].re = ac * fn23 - pp * tfn_im;
				up[1].im = ac * tfn_im + pp * fn23;
				bsum->re = up[1].re + zc.re;
				bsum->im = up[1].im + zc.im;

				if(rfnu >= tol)
				{
					zth_re = rzth_re;
					zth_im = rzth_im;
					w2.re = fn23;
					w2.im = tfn_im;
					up[0].re = 1.0;
					up[0].im = 0.0;
					pp = 1.0;
					btol = tol * (fabs(bsum->re) + fabs(bsum->im));
					ks = -1;
					kp1 = 2;
					l = 2;
					ias = 0;
					ibs = 0;

					for(l1 = 0; l1 < 14; l1++)
					{
						cr[l1].re = 0.0;
						cr[l1].im = 0.0;
						dr[l1].re = 0.0;
						dr[l1].im = 0.0;
					}

					lr = 2;
					exitg5 = false;

					while((!exitg5) && (lr < 13))
					{
						for(k = lr; k <= lr + 1; k++)
						{
							ks++;
							kp1++;
							l++;
							suma.re = C[l];
							suma.im = 0.0;

							for(l1 = 2; l1 <= kp1; l1++)
							{
								l++;
								ac = suma.re;
								suma.re = (suma.re * w.re - suma.im * w.im) + C[l];
								suma.im = ac * w.im + suma.im * w.re;
							}

							ac = w2.re;
							w2.re = w2.re * fn23 - w2.im * tfn_im;
							w2.im = ac * tfn_im + w2.im * fn23;
							up[kp1 - 1].re = w2.re * suma.re - w2.im * suma.im;
							up[kp1 - 1].im = w2.re * suma.im + w2.im * suma.re;
							cr[ks].re = BR[ks + 1] * zth_re;
							cr[ks].im = BR[ks + 1] * zth_im;
							ac = zth_re;
							zth_re = zth_re * rzth_re - zth_im * rzth_im;
							zth_im = ac * rzth_im + zth_im * rzth_re;
							dr[ks].re = AR[ks + 2] * zth_re;
							dr[ks].im = AR[ks + 2] * zth_im;
						}

						pp *= rfnu2;

						if(ias != 1)
						{
							suma = up[lr];
							l1 = lr;

							for(l2 = 1; l2 <= lr; l2++)
							{
								l1--;
								suma.re += cr[l2 - 1].re * up[l1].re - cr[l2 - 1].im * up[l1].im;
								suma.im += cr[l2 - 1].re * up[l1].im + cr[l2 - 1].im * up[l1].re;
							}

							asum->re += suma.re;
							asum->im += suma.im;

							if((pp < tol) && (fabs(asum->re) + fabs(asum->im) < tol))
								ias = 1;
						}

						if(ibs != 1)
						{
							zb.re = up[lr + 1].re + (up[lr].re * zc.re - up[lr].im * zc.im);
							zb.im = up[lr + 1].im + (up[lr].re * zc.im + up[lr].im * zc.re);
							l1 = lr;

							for(l2 = 1; l2 <= lr; l2++)
							{
								l1--;
								zb.re += dr[l2 - 1].re * up[l1].re - dr[l2 - 1].im * up[l1].im;
								zb.im += dr[l2 - 1].re * up[l1].im + dr[l2 - 1].im * up[l1].re;
							}

							bsum->re += zb.re;
							bsum->im += zb.im;

							if((pp < btol) && (fabs(bsum->re) + fabs(bsum->im) < tol))
								ibs = 1;
						}

						if((ias == 1) && (ibs == 1))
							exitg5 = true;
						else
							lr += 2;
					}
				}

				asum->re++;
				bsum->re = -bsum->re;
				bsum->im = -bsum->im;
				pp = bsum->re * rfn13_re - bsum->im * 0.0;
				ac = bsum->re * 0.0 + bsum->im * rfn13_re;

				if(rtzta_im == 0.0)
				{
					if(ac == 0.0)
					{
						bsum->re = pp / rtzta_re;
						bsum->im = 0.0;
					}
					else if(pp == 0.0)
					{
						bsum->re = 0.0;
						bsum->im = ac / rtzta_re;
					}
					else
					{
						bsum->re = pp / rtzta_re;
						bsum->im = ac / rtzta_re;
					}
				}
				else if(rtzta_re == 0.0)
				{
					if(pp == 0.0)
					{
						bsum->re = ac / rtzta_im;
						bsum->im = 0.0;
					}
					else if(ac == 0.0)
					{
						bsum->re = 0.0;
						bsum->im = -(pp / rtzta_im);
					}
					else
					{
						bsum->re = ac / rtzta_im;
						bsum->im = -(pp / rtzta_im);
					}
				}
				else
				{
					brm = fabs(rtzta_re);
					btol = fabs(rtzta_im);

					if(brm > btol)
					{
						btol = rtzta_im / rtzta_re;
						d = rtzta_re + btol * rtzta_im;
						bsum->re = (pp + btol * ac) / d;
						bsum->im = (ac - btol * pp) / d;
					}
					else if(btol == brm)
					{
						if(rtzta_re > 0.0)
							btol = 0.5;
						else
							btol = -0.5;

						if(rtzta_im > 0.0)
							d = 0.5;
						else
							d = -0.5;

						bsum->re = (pp * btol + ac * d) / brm;
						bsum->im = (ac * btol - pp * d) / brm;
					}
					else
					{
						btol = rtzta_re / rtzta_im;
						d = rtzta_im + btol * rtzta_re;
						bsum->re = (btol * pp + ac) / d;
						bsum->im = (btol * ac - pp) / d;
					}
				}
			}
		}
		else
		{
			k = 1;
			memset(&p[0], 0, 30U * sizeof(complex));
			memset(&ap[0], 0, 30U * sizeof(double));
			p[0].re = 1.0;
			p[0].im = 0.0;
			suma.re = 0.6299605249474366;
			suma.im = 0.0;
			ap[0] = 1.0;

			if(ac >= tol)
			{
				k = 30;
				l1 = 1;
				exitg4 = false;

				while((!exitg4) && (l1 + 1 < 31))
				{
					pp = p[l1 - 1].re;
					btol = p[l1 - 1].im;
					p[l1].re = pp * w2.re - btol * w2.im;
					p[l1].im = pp * w2.im + btol * w2.re;
					suma.re += p[l1].re * GAMA[l1];
					suma.im += p[l1].im * GAMA[l1];
					ap[l1] = ap[l1 - 1] * ac;

					if(ap[l1] < tol)
					{
						k = l1 + 1;
						exitg4 = true;
					}
					else
						l1++;
				}
			}

			zb.re = w2.re * suma.re - w2.im * suma.im;
			zb.im = w2.re * suma.im + w2.im * suma.re;
			arg->re = fn23 * zb.re;
			arg->im = fn23 * zb.im;
			b_sqrt(&suma);
			b_sqrt(&w2);
			zeta2->re = fnu * w2.re;
			zeta2->im = fnu * w2.im;
			ac = (zb.re * suma.re - zb.im * suma.im) * 0.66666666666666663 + 1.0;
			pp = (zb.re * suma.im + zb.im * suma.re) * 0.66666666666666663;
			zeta1->re = zeta2->re * ac - zeta2->im * pp;
			zeta1->im = zeta2->re * pp + zeta2->im * ac;
			suma.re += suma.re;
			suma.im += suma.im;
			b_sqrt(&suma);
			phi->re = suma.re * rfn13_re - suma.im * 0.0;
			phi->im = suma.re * 0.0 + suma.im * rfn13_re;

			if(ipmtr == 1)
				;
			else
			{
				zb.re = 0.0;
				zb.im = 0.0;

				for(l1 = 0; l1 + 1 <= k; l1++)
				{
					zb.re += p[l1].re * BETA[l1];
					zb.im += p[l1].im * BETA[l1];
				}

				*bsum = zb;
				l1 = 0;
				l2 = 30;
				btol = tol * rt_hypotd_snf(zb.re, zb.im);
				ac = tol;
				pp = 1.0;
				ias = 0;
				ibs = 0;

				if(rfnu2 >= tol)
				{
					lr = 2;
					exitg1 = false;

					while((!exitg1) && (lr < 8))
					{
						ac /= rfnu2;
						pp *= rfnu2;

						if(ias != 1)
						{
							suma.re = 0.0;
							suma.im = 0.0;
							ks = 0;
							exitg3 = false;

							while((!exitg3) && (ks + 1 <= k))
							{
								suma.re += p[ks].re * ALFA[l1 + ks];
								suma.im += p[ks].im * ALFA[l1 + ks];

								if(ap[ks] < ac)
									exitg3 = true;
								else
									ks++;
							}

							asum->re += suma.re * pp;
							asum->im += suma.im * pp;

							if(pp < tol)
								ias = 1;
						}

						if(ibs != 1)
						{
							zb.re = 0.0;
							zb.im = 0.0;
							ks = 0;
							exitg2 = false;

							while((!exitg2) && (ks + 1 <= k))
							{
								zb.re += p[ks].re * BETA[l2 + ks];
								zb.im += p[ks].im * BETA[l2 + ks];

								if(ap[ks] < ac)
									exitg2 = true;
								else
									ks++;
							}

							bsum->re += zb.re * pp;
							bsum->im += zb.im * pp;

							if(pp < btol)
								ibs = 1;
						}

						if((ias == 1) && (ibs == 1))
							exitg1 = true;
						else
						{
							l1 += 30;
							l2 += 30;
							lr++;
						}
					}
				}

				asum->re++;
				bsum->re *= rfnu * rfn13_re;
				bsum->im *= rfnu * rfn13_re;
			}
		}
	}
}

void cunhj(complex z, double fnu, double tol, complex *phi, complex *arg, complex *zeta1, complex *zeta2)
{
	double ac;
	complex zb;
	double fn23;
	double rfn13_re;
	complex w2;
	complex p[30];
	double ap[30];
	double tmpr;
	double tmpi;
	complex suma;
	int i;
	bool exitg1;
	double brm;
	complex zc;
	double p_im;
	static const double GAMA[30] = { 0.6299605249474366, 0.25198420997897464,
		0.15479030041565583, 0.11071306241615901, 0.085730939552739485,
		0.069716131695868433, 0.058608567189371359, 0.050469887353631067,
		0.044260058068915482, 0.039372066154350994, 0.035428319592445537,
		0.032181885750209825, 0.029464624079115768, 0.027158167711293448,
		0.025176827297386177, 0.023457075530607888, 0.021950839013490719,
		0.020621082823564625, 0.019438824089788084, 0.018381063380068317,
		0.017429321323196318, 0.016568583778661234, 0.015786528598791844,
		0.015072950149409559, 0.014419325083995464, 0.013818480573534178,
		0.013264337899427657, 0.012751712197049864, 0.012276154531876277,
		0.011833826239848241 };

	double zth_re;
	double zth_im;
	ac = fnu * 2.2250738585072014E-305;

	if((fabs(z.re) <= ac) && (fabs(z.im) <= ac))
	{
		zeta1->re = 1402.9773265065639 + fnu;
		zeta1->im = 0.0;
		zeta2->re = fnu;
		zeta2->im = 0.0;
		phi->re = 1.0;
		phi->im = 0.0;
		arg->re = 1.0;
		arg->im = 0.0;
	}
	else
	{
		zb.re = 1.0 / fnu * z.re;
		zb.im = 1.0 / fnu * z.im;
		ac = rt_powd_snf(fnu, 0.33333333333333331);
		fn23 = ac * ac;
		rfn13_re = 1.0 / ac;
		w2.re = 1.0 - (zb.re * zb.re - zb.im * zb.im);
		w2.im = 0.0 - (zb.re * zb.im + zb.im * zb.re);
		ac = rt_hypotd_snf(w2.re, w2.im);

		if(ac > 0.25)
		{
			b_sqrt(&w2);
			tmpr = w2.re;
			tmpi = w2.im;

			if(w2.re < 0.0)
				tmpr = 0.0;

			if(w2.im < 0.0)
				tmpi = 0.0;

			w2.re = tmpr;
			w2.im = tmpi;

			if(zb.im == 0.0)
			{
				if(tmpi == 0.0)
				{
					zc.re = (1.0 + tmpr) / zb.re;
					zc.im = 0.0;
				}
				else if(1.0 + tmpr == 0.0)
				{
					zc.re = 0.0;
					zc.im = tmpi / zb.re;
				}
				else
				{
					zc.re = (1.0 + tmpr) / zb.re;
					zc.im = tmpi / zb.re;
				}
			}
			else if(zb.re == 0.0)
			{
				if(1.0 + tmpr == 0.0)
				{
					zc.re = tmpi / zb.im;
					zc.im = 0.0;
				}
				else if(tmpi == 0.0)
				{
					zc.re = 0.0;
					zc.im = -((1.0 + tmpr) / zb.im);
				}
				else
				{
					zc.re = tmpi / zb.im;
					zc.im = -((1.0 + tmpr) / zb.im);
				}
			}
			else
			{
				brm = fabs(zb.re);
				ac = fabs(zb.im);

				if(brm > ac)
				{
					ac = zb.im / zb.re;
					p_im = zb.re + ac * zb.im;
					zc.re = ((1.0 + tmpr) + ac * tmpi) / p_im;
					zc.im = (tmpi - ac * (1.0 + tmpr)) / p_im;
				}
				else if(ac == brm)
				{
					if(zb.re > 0.0)
						ac = 0.5;
					else
						ac = -0.5;

					if(zb.im > 0.0)
						p_im = 0.5;
					else
						p_im = -0.5;

					zc.re = ((1.0 + tmpr) * ac + tmpi * p_im) / brm;
					zc.im = (tmpi * ac - (1.0 + tmpr) * p_im) / brm;
				}
				else
				{
					ac = zb.re / zb.im;
					p_im = zb.im + ac * zb.re;
					zc.re = (ac * (1.0 + tmpr) + tmpi) / p_im;
					zc.im = (ac * tmpi - (1.0 + tmpr)) / p_im;
				}
			}

			b_log(&zc);
			tmpr = zc.re;
			tmpi = zc.im;

			if(zc.re < 0.0)
				tmpr = 0.0;

			if(zc.im < 0.0)
				tmpi = 0.0;
			else
			{
				if(zc.im > 1.5707963267948966)
					tmpi = 1.5707963267948966;
			}

			zth_re = 1.5 * (tmpr - w2.re);
			zth_im = 1.5 * (tmpi - w2.im);
			zeta1->re = fnu * tmpr;
			zeta1->im = fnu * tmpi;
			zeta2->re = fnu * w2.re;
			zeta2->im = fnu * w2.im;

			if((zth_re >= 0.0) && (zth_im < 0.0))
				ac = 4.71238898038469;
			else if(zth_re != 0.0)
			{
				ac = atan(zth_im / zth_re);

				if(zth_re < 0.0)
					ac += 3.1415926535897931;
			}
			else
				ac = 1.5707963267948966;

			tmpr = rt_powd_snf(rt_hypotd_snf(zth_re, zth_im), 0.66666666666666663);
			ac *= 0.66666666666666663;
			zb.re = tmpr * cos(ac);
			zb.im = tmpr * sin(ac);

			if(zb.im < 0.0)
				zb.im = 0.0;

			arg->re = fn23 * zb.re;
			arg->im = fn23 * zb.im;

			if(zb.im == 0.0)
			{
				if(zth_im == 0.0)
				{
					tmpr = zth_re / zb.re;
					zth_im = 0.0;
				}
				else if(zth_re == 0.0)
				{
					tmpr = 0.0;
					zth_im /= zb.re;
				}
				else
				{
					tmpr = zth_re / zb.re;
					zth_im /= zb.re;
				}
			}
			else if(zb.re == 0.0)
			{
				if(zth_re == 0.0)
				{
					tmpr = zth_im / zb.im;
					zth_im = 0.0;
				}
				else if(zth_im == 0.0)
				{
					tmpr = 0.0;
					zth_im = -(zth_re / zb.im);
				}
				else
				{
					tmpr = zth_im / zb.im;
					zth_im = -(zth_re / zb.im);
				}
			}
			else
			{
				brm = fabs(zb.re);
				ac = zb.im;

				if(brm > ac)
				{
					ac = zb.im / zb.re;
					p_im = zb.re + ac * zb.im;
					tmpr = (zth_re + ac * zth_im) / p_im;
					zth_im = (zth_im - ac * zth_re) / p_im;
				}
				else if(ac == brm)
				{
					if(zb.re > 0.0)
						ac = 0.5;
					else
						ac = -0.5;

					if(zb.im > 0.0)
						p_im = 0.5;
					else
						p_im = -0.5;

					tmpr = (zth_re * ac + zth_im * p_im) / brm;
					zth_im = (zth_im * ac - zth_re * p_im) / brm;
				}
				else
				{
					ac = zb.re / zb.im;
					p_im = zb.im + ac * zb.re;
					tmpr = (ac * zth_re + zth_im) / p_im;
					zth_im = (ac * zth_im - zth_re) / p_im;
				}
			}

			if(w2.im == 0.0)
			{
				if(zth_im == 0.0)
				{
					suma.re = tmpr / w2.re;
					suma.im = 0.0;
				}
				else if(tmpr == 0.0)
				{
					suma.re = 0.0;
					suma.im = zth_im / w2.re;
				}
				else
				{
					suma.re = tmpr / w2.re;
					suma.im = zth_im / w2.re;
				}
			}
			else if(w2.re == 0.0)
			{
				if(tmpr == 0.0)
				{
					suma.re = zth_im / w2.im;
					suma.im = 0.0;
				}
				else if(zth_im == 0.0)
				{
					suma.re = 0.0;
					suma.im = -(tmpr / w2.im);
				}
				else
				{
					suma.re = zth_im / w2.im;
					suma.im = -(tmpr / w2.im);
				}
			}
			else
			{
				brm = fabs(w2.re);
				ac = fabs(w2.im);

				if(brm > ac)
				{
					ac = w2.im / w2.re;
					p_im = w2.re + ac * w2.im;
					suma.re = (tmpr + ac * zth_im) / p_im;
					suma.im = (zth_im - ac * tmpr) / p_im;
				}
				else if(ac == brm)
				{
					if(w2.re > 0.0)
						ac = 0.5;
					else
						ac = -0.5;

					if(w2.im > 0.0)
						p_im = 0.5;
					else
						p_im = -0.5;

					suma.re = (tmpr * ac + zth_im * p_im) / brm;
					suma.im = (zth_im * ac - tmpr * p_im) / brm;
				}
				else
				{
					ac = w2.re / w2.im;
					p_im = w2.im + ac * w2.re;
					suma.re = (ac * tmpr + zth_im) / p_im;
					suma.im = (ac * zth_im - tmpr) / p_im;
				}
			}

			zb.re = suma.re + suma.re;
			zb.im = suma.im + suma.im;
			b_sqrt(&zb);
			phi->re = zb.re * rfn13_re - zb.im * 0.0;
			phi->im = zb.re * 0.0 + zb.im * rfn13_re;
		}
		else
		{
			memset(&p[0], 0, 30U * sizeof(complex));
			memset(&ap[0], 0, 30U * sizeof(double));
			p[0].re = 1.0;
			p[0].im = 0.0;
			suma.re = 0.6299605249474366;
			suma.im = 0.0;
			ap[0] = 1.0;

			if(ac >= tol)
			{
				i = 1;
				exitg1 = false;

				while((!exitg1) && (i + 1 < 31))
				{
					tmpr = p[i - 1].re;
					p_im = p[i - 1].im;
					p[i].re = tmpr * w2.re - p_im * w2.im;
					p[i].im = tmpr * w2.im + p_im * w2.re;
					suma.re += p[i].re * GAMA[i];
					suma.im += p[i].im * GAMA[i];
					ap[i] = ap[i - 1] * ac;

					if(ap[i] < tol)
						exitg1 = true;
					else
						i++;
				}
			}

			zb.re = w2.re * suma.re - w2.im * suma.im;
			zb.im = w2.re * suma.im + w2.im * suma.re;
			arg->re = fn23 * zb.re;
			arg->im = fn23 * zb.im;
			b_sqrt(&suma);
			b_sqrt(&w2);
			zeta2->re = fnu * w2.re;
			zeta2->im = fnu * w2.im;
			ac = (zb.re * suma.re - zb.im * suma.im) * 0.66666666666666663 + 1.0;
			tmpr = (zb.re * suma.im + zb.im * suma.re) * 0.66666666666666663;
			zeta1->re = zeta2->re * ac - zeta2->im * tmpr;
			zeta1->im = zeta2->re * tmpr + zeta2->im * ac;
			suma.re += suma.re;
			suma.im += suma.im;
			b_sqrt(&suma);
			phi->re = suma.re * rfn13_re - suma.im * 0.0;
			phi->im = suma.re * 0.0 + suma.im * rfn13_re;
		}
	}
}

void b_cunik(complex zr, double fnu, int ikflg, int ipmtr, double tol, int *init, complex cwrk[16], complex *phi, complex *zeta1, complex *zeta2, complex *summ)
{
	double rfn;
	double ac;
	bool guard1 = false;
	double t_re;
	double t_im;
	double s_re;
	double s_im;
	complex sr;
	int i;
	double b_t_re;
	double cfn_re;
	double ai;
	double crfn_re;
	int l;
	int exitg1;
	int j;
	static const double C[120] = { 1.0, -0.20833333333333334, 0.125,
		0.3342013888888889, -0.40104166666666669, 0.0703125, -1.0258125964506173,
		1.8464626736111112, -0.8912109375, 0.0732421875, 4.6695844234262474,
		-11.207002616222994, 8.78912353515625, -2.3640869140625, 0.112152099609375,
		-28.212072558200244, 84.636217674600729, -91.818241543240021,
		42.534998745388457, -7.3687943594796321, 0.22710800170898438,
		212.57013003921713, -765.25246814118168, 1059.9904525279999,
		-699.57962737613252, 218.19051174421159, -26.491430486951554,
		0.57250142097473145, -1919.4576623184071, 8061.7221817373093,
		-13586.550006434138, 11655.393336864534, -5305.646978613403,
		1200.9029132163525, -108.09091978839466, 1.7277275025844574,
		20204.291330966149, -96980.598388637518, 192547.00123253153,
		-203400.17728041555, 122200.46498301746, -41192.65496889755,
		7109.5143024893641, -493.915304773088, 6.074042001273483,
		-242919.18790055133, 1.3117636146629772E+6, -2.9980159185381066E+6,
		3.7632712976564039E+6, -2.8135632265865342E+6, 1.2683652733216248E+6,
		-331645.17248456361, 45218.768981362729, -2499.8304818112097,
		24.380529699556064, 3.2844698530720379E+6, -1.9706819118432228E+7,
		5.0952602492664643E+7, -7.4105148211532652E+7, 6.6344512274729028E+7,
		-3.7567176660763353E+7, 1.3288767166421818E+7, -2.7856181280864547E+6,
		308186.40461266239, -13886.08975371704, 110.01714026924674,
		-4.932925366450996E+7, 3.2557307418576574E+8, -9.394623596815784E+8,
		1.55359689957058E+9, -1.6210805521083372E+9, 1.1068428168230145E+9,
		-4.9588978427503031E+8, 1.4206290779753309E+8, -2.447406272573873E+7,
		2.2437681779224495E+6, -84005.433603024081, 551.33589612202059,
		8.1478909611831212E+8, -5.8664814920518475E+9, 1.8688207509295826E+10,
		-3.4632043388158775E+10, 4.1280185579753975E+10, -3.3026599749800724E+10,
		1.79542137311556E+10, -6.5632937926192846E+9, 1.5592798648792574E+9,
		-2.2510566188941526E+8, 1.7395107553978164E+7, -549842.32757228869,
		3038.0905109223841, -1.4679261247695616E+10, 1.144982377320258E+11,
		-3.9909617522446649E+11, 8.1921866954857727E+11, -1.0983751560812233E+12,
		1.0081581068653821E+12, -6.4536486924537646E+11, 2.8790064990615057E+11,
		-8.786707217802327E+10, 1.7634730606834969E+10, -2.1671649832237949E+9,
		1.4315787671888897E+8, -3.8718334425726128E+6, 18257.755474293175,
		2.86464035717679E+11, -2.4062979000285039E+12, 9.1093411852398984E+12,
		-2.0516899410934437E+13, 3.056512551993532E+13, -3.166708858478516E+13,
		2.334836404458184E+13, -1.2320491305598287E+13, 4.6127257808491318E+12,
		-1.1965528801961816E+12, 2.0591450323241E+11, -2.1822927757529224E+10,
		1.2470092935127103E+9, -2.9188388122220814E+7, 118838.42625678325 };

	rfn = 1.0 / fnu;
	ac = fnu * 2.2250738585072014E-305;
	guard1 = false;

	if(*init == 0)
	{
		summ->re = 0.0;
		summ->im = 0.0;

		if((fabs(zr.re) > ac) || (fabs(zr.im) > ac))
		{
			t_re = zr.re * rfn - zr.im * 0.0;
			t_im = zr.re * 0.0 + zr.im * rfn;
			s_re = (t_re * t_re - t_im * t_im) + 1.0;
			s_im = t_re * t_im + t_im * t_re;
			sr.re = s_re;
			sr.im = s_im;
			b_sqrt(&sr);
			b_t_re = t_re;
			cfn_re = 1.0 + sr.re;
			ai = sr.im;

			if(t_im == 0.0)
			{
				if(ai == 0.0)
				{
					t_re = cfn_re / t_re;
					t_im = 0.0;
				}
				else if(cfn_re == 0.0)
				{
					t_re = 0.0;
					t_im = ai / b_t_re;
				}
				else
				{
					t_re = cfn_re / t_re;
					t_im = ai / b_t_re;
				}
			}
			else if(t_re == 0.0)
			{
				if(cfn_re == 0.0)
				{
					t_re = ai / t_im;
					t_im = 0.0;
				}
				else if(ai == 0.0)
				{
					t_re = 0.0;
					t_im = -(cfn_re / t_im);
				}
				else
				{
					t_re = ai / t_im;
					t_im = -(cfn_re / t_im);
				}
			}
			else
			{
				crfn_re = fabs(t_re);
				b_t_re = fabs(t_im);

				if(crfn_re > b_t_re)
				{
					ac = t_im / t_re;
					b_t_re = t_re + ac * t_im;
					t_re = (cfn_re + ac * ai) / b_t_re;
					t_im = (ai - ac * cfn_re) / b_t_re;
				}
				else if (b_t_re == crfn_re)
				{
					if(t_re > 0.0)
						ac = 0.5;
					else
						ac = -0.5;

					if(t_im > 0.0)
						b_t_re = 0.5;
					else
						b_t_re = -0.5;

					t_re = (cfn_re * ac + ai * b_t_re) / crfn_re;
					t_im = (ai * ac - cfn_re * b_t_re) / crfn_re;
				}
				else
				{
					ac = t_re / t_im;
					b_t_re = t_im + ac * t_re;
					t_re = (ac * cfn_re + ai) / b_t_re;
					t_im = (ac * ai - cfn_re) / b_t_re;
				}
			}

			if((t_im == 0.0) && rtIsNaN(t_re))
				;
			else if((fabs(t_re) > 8.9884656743115785E+307) || (fabs(t_im) > 8.9884656743115785E+307))
			{
				b_t_re = t_re;
				t_re = log(rt_hypotd_snf(t_re / 2.0, t_im / 2.0)) + 0.69314718055994529;
				t_im = rt_atan2d_snf(t_im, b_t_re);
			}
			else
			{
				b_t_re = t_re;
				t_re = log(rt_hypotd_snf(t_re, t_im));
				t_im = rt_atan2d_snf(t_im, b_t_re);
			}

			zeta1->re = fnu * t_re - 0.0 * t_im;
			zeta1->im = fnu * t_im + 0.0 * t_re;
			zeta2->re = fnu * sr.re - 0.0 * sr.im;
			zeta2->im = fnu * sr.im + 0.0 * sr.re;

			if(sr.im == 0.0)
			{
				sr.re = rfn / sr.re;
				sr.im = 0.0;
			}
			else if (sr.re == 0.0)
			{
				if(rfn == 0.0)
				{
					sr.re = 0.0 / sr.im;
					sr.im = 0.0;
				}
				else
				{
					sr.re = 0.0;
					sr.im = -(rfn / sr.im);
				}
			}
			else
			{
				crfn_re = fabs(sr.re);
				b_t_re = fabs(sr.im);

				if(crfn_re > b_t_re)
				{
					ac = sr.im / sr.re;
					b_t_re = sr.re + ac * sr.im;
					sr.re = (rfn + ac * 0.0) / b_t_re;
					sr.im = (0.0 - ac * rfn) / b_t_re;
				}
				else if(b_t_re == crfn_re)
				{
					if(sr.re > 0.0)
						ac = 0.5;
					else
						ac = -0.5;

					if(sr.im > 0.0)
						b_t_re = 0.5;
					else
						b_t_re = -0.5;

					sr.re = (rfn * ac + 0.0 * b_t_re) / crfn_re;
					sr.im = (0.0 * ac - rfn * b_t_re) / crfn_re;
				}
				else
				{
					ac = sr.re / sr.im;
					b_t_re = sr.im + ac * sr.re;
					sr.re = ac * rfn / b_t_re;
					sr.im = (ac * 0.0 - rfn) / b_t_re;
				}
			}

			cwrk[15] = sr;
			b_sqrt(&cwrk[15]);
			phi->re = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) * cwrk[15].re;
			phi->im = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) * cwrk[15].im;

			if(ipmtr != 0)
				;
			else
			{
				if(s_im == 0.0)
				{
					cfn_re = 1.0 / s_re;
					t_im = 0.0;
				}
				else if(s_re == 0.0)
				{
					cfn_re = 0.0;
					t_im = -(1.0 / s_im);
				}
				else
				{
					crfn_re = fabs(s_re);
					b_t_re = fabs(s_im);

					if(crfn_re > b_t_re)
					{
						ac = s_im / s_re;
						b_t_re = s_re + ac * s_im;
						cfn_re = (1.0 + ac * 0.0) / b_t_re;
						t_im = (0.0 - ac) / b_t_re;
					}
					else if(b_t_re == crfn_re)
					{
						if(s_re > 0.0)
							ac = 0.5;
						else
							ac = -0.5;

						if(s_im > 0.0)
							b_t_re = 0.5;
						else
							b_t_re = -0.5;

						cfn_re = (ac + 0.0 * b_t_re) / crfn_re;
						t_im = (0.0 * ac - b_t_re) / crfn_re;
					}
					else
					{
						ac = s_re / s_im;
						b_t_re = s_im + ac * s_re;
						cfn_re = ac / b_t_re;
						t_im = (ac * 0.0 - 1.0) / b_t_re;
					}
				}

				cwrk[0].re = 1.0;
				cwrk[0].im = 0.0;
				crfn_re = 1.0;
				t_re = 0.0;
				ac = 1.0;
				l = 0;
				*init = 15;
				i = 1;

				do
				{
					exitg1 = 0;

					if(i < 15)
					{
						s_re = 0.0;
						s_im = 0.0;

						for(j = 0; j <= i; j++)
						{
							l++;
							b_t_re = s_re;
							s_re = (s_re * cfn_re - s_im * t_im) + C[l];
							s_im = b_t_re * t_im + s_im * cfn_re;
						}

						b_t_re = crfn_re;
						crfn_re = crfn_re * sr.re - t_re * sr.im;
						t_re = b_t_re * sr.im + t_re * sr.re;
						cwrk[i].re = crfn_re * s_re - t_re * s_im;
						cwrk[i].im = crfn_re * s_im + t_re * s_re;
						ac *= rfn;

						if((ac < tol) && (fabs(cwrk[i].re) + fabs(cwrk[i].im) < tol))
						{
							*init = i + 1;
							guard1 = true;
							exitg1 = 1;
						}
						else
							i++;
					}
					else
					{
						guard1 = true;
						exitg1 = 1;
					}
				}
				while(exitg1 == 0);
			}
		}
		else
		{
			zeta1->re = 1402.9773265065639 + fnu;
			zeta1->im = 0.0;
			zeta2->re = fnu;
			zeta2->im = 0.0;
			phi->re = 1.0;
			phi->im = 0.0;
		}
	}
	else
	{
		zeta1->re = 0.0;
		zeta1->im = 0.0;
		zeta2->re = 0.0;
		zeta2->im = 0.0;
		guard1 = true;
	}

	if(guard1)
	{
		if(ikflg == 2)
		{
			s_re = 0.0;
			s_im = 0.0;
			t_re = 1.0;
			t_im = 0.0;

			for(i = 1; i <= *init; i++)
			{
				s_re += t_re * cwrk[i - 1].re - t_im * cwrk[i - 1].im;
				s_im += t_re * cwrk[i - 1].im + t_im * cwrk[i - 1].re;
				t_re = -t_re;
				t_im = -t_im;
			}

			summ->re = s_re;
			summ->im = s_im;
			phi->re = 1.2533141373155003 * cwrk[15].re;
			phi->im = 1.2533141373155003 * cwrk[15].im;
		}
		else
		{
			s_re = 0.0;
			s_im = 0.0;

			for(i = 1; i <= *init; i++)
			{
				s_re += cwrk[i - 1].re;
				s_im += cwrk[i - 1].im;
			}

			summ->re = s_re;
			summ->im = s_im;
			phi->re = 0.3989422804014327 * cwrk[15].re;
			phi->im = 0.3989422804014327 * cwrk[15].im;
		}
	}
}

void cunik(complex zr, double fnu, int ikflg, int ipmtr, double tol, int init, complex *phi, complex *zeta1, complex *zeta2)
{
	complex cwrk[16];
	int i;
	double rfn;
	double ac;
	bool guard1 = false;
	double t_re;
	double t_im;
	double s_re;
	double s_im;
	complex sr;
	double b_t_re;
	double cfn_re;
	double ai;
	double crfn_re;
	int l;
	int exitg1;
	int j;
	static const double C[120] = { 1.0, -0.20833333333333334, 0.125,
		0.3342013888888889, -0.40104166666666669, 0.0703125, -1.0258125964506173,
		1.8464626736111112, -0.8912109375, 0.0732421875, 4.6695844234262474,
		-11.207002616222994, 8.78912353515625, -2.3640869140625, 0.112152099609375,
		-28.212072558200244, 84.636217674600729, -91.818241543240021,
		42.534998745388457, -7.3687943594796321, 0.22710800170898438,
		212.57013003921713, -765.25246814118168, 1059.9904525279999,
		-699.57962737613252, 218.19051174421159, -26.491430486951554,
		0.57250142097473145, -1919.4576623184071, 8061.7221817373093,
		-13586.550006434138, 11655.393336864534, -5305.646978613403,
		1200.9029132163525, -108.09091978839466, 1.7277275025844574,
		20204.291330966149, -96980.598388637518, 192547.00123253153,
		-203400.17728041555, 122200.46498301746, -41192.65496889755,
		7109.5143024893641, -493.915304773088, 6.074042001273483,
		-242919.18790055133, 1.3117636146629772E+6, -2.9980159185381066E+6,
		3.7632712976564039E+6, -2.8135632265865342E+6, 1.2683652733216248E+6,
		-331645.17248456361, 45218.768981362729, -2499.8304818112097,
		24.380529699556064, 3.2844698530720379E+6, -1.9706819118432228E+7,
		5.0952602492664643E+7, -7.4105148211532652E+7, 6.6344512274729028E+7,
		-3.7567176660763353E+7, 1.3288767166421818E+7, -2.7856181280864547E+6,
		308186.40461266239, -13886.08975371704, 110.01714026924674,
		-4.932925366450996E+7, 3.2557307418576574E+8, -9.394623596815784E+8,
		1.55359689957058E+9, -1.6210805521083372E+9, 1.1068428168230145E+9,
		-4.9588978427503031E+8, 1.4206290779753309E+8, -2.447406272573873E+7,
		2.2437681779224495E+6, -84005.433603024081, 551.33589612202059,
		8.1478909611831212E+8, -5.8664814920518475E+9, 1.8688207509295826E+10,
		-3.4632043388158775E+10, 4.1280185579753975E+10, -3.3026599749800724E+10,
		1.79542137311556E+10, -6.5632937926192846E+9, 1.5592798648792574E+9,
		-2.2510566188941526E+8, 1.7395107553978164E+7, -549842.32757228869,
		3038.0905109223841, -1.4679261247695616E+10, 1.144982377320258E+11,
		-3.9909617522446649E+11, 8.1921866954857727E+11, -1.0983751560812233E+12,
		1.0081581068653821E+12, -6.4536486924537646E+11, 2.8790064990615057E+11,
		-8.786707217802327E+10, 1.7634730606834969E+10, -2.1671649832237949E+9,
		1.4315787671888897E+8, -3.8718334425726128E+6, 18257.755474293175,
		2.86464035717679E+11, -2.4062979000285039E+12, 9.1093411852398984E+12,
		-2.0516899410934437E+13, 3.056512551993532E+13, -3.166708858478516E+13,
		2.334836404458184E+13, -1.2320491305598287E+13, 4.6127257808491318E+12,
		-1.1965528801961816E+12, 2.0591450323241E+11, -2.1822927757529224E+10,
		1.2470092935127103E+9, -2.9188388122220814E+7, 118838.42625678325 };

	for(i = 0; i < 16; i++)
	{
		cwrk[i].re = 0.0;
		cwrk[i].im = 0.0;
	}

	rfn = 1.0 / fnu;
	ac = fnu * 2.2250738585072014E-305;
	guard1 = false;

	if(init == 0)
	{
		if((fabs(zr.re) > ac) || (fabs(zr.im) > ac))
		{
			t_re = zr.re * rfn - zr.im * 0.0;
			t_im = zr.re * 0.0 + zr.im * rfn;
			s_re = (t_re * t_re - t_im * t_im) + 1.0;
			s_im = t_re * t_im + t_im * t_re;
			sr.re = s_re;
			sr.im = s_im;
			b_sqrt(&sr);
			b_t_re = t_re;
			cfn_re = 1.0 + sr.re;
			ai = sr.im;

			if(t_im == 0.0)
			{
				if(ai == 0.0)
				{
					t_re = cfn_re / t_re;
					t_im = 0.0;
				}
				else if (cfn_re == 0.0)
				{
					t_re = 0.0;
					t_im = ai / b_t_re;
				}
				else
				{
					t_re = cfn_re / t_re;
					t_im = ai / b_t_re;
				}
			}
			else if(t_re == 0.0)
			{
				if(cfn_re == 0.0)
				{
					t_re = ai / t_im;
					t_im = 0.0;
				}
				else if(ai == 0.0)
				{
					t_re = 0.0;
					t_im = -(cfn_re / t_im);
				}
				else
				{
					t_re = ai / t_im;
					t_im = -(cfn_re / t_im);
				}
			}
			else
			{
				crfn_re = fabs(t_re);
				b_t_re = fabs(t_im);

				if(crfn_re > b_t_re)
				{
					ac = t_im / t_re;
					b_t_re = t_re + ac * t_im;
					t_re = (cfn_re + ac * ai) / b_t_re;
					t_im = (ai - ac * cfn_re) / b_t_re;
				}
				else if(b_t_re == crfn_re)
				{
					if(t_re > 0.0)
						ac = 0.5;
					else
						ac = -0.5;

					if(t_im > 0.0)
						b_t_re = 0.5;
					else
						b_t_re = -0.5;

					t_re = (cfn_re * ac + ai * b_t_re) / crfn_re;
					t_im = (ai * ac - cfn_re * b_t_re) / crfn_re;
				}
				else
				{
					ac = t_re / t_im;
					b_t_re = t_im + ac * t_re;
					t_re = (ac * cfn_re + ai) / b_t_re;
					t_im = (ac * ai - cfn_re) / b_t_re;
				}
			}

			if((t_im == 0.0) && rtIsNaN(t_re))
				;
			else if((fabs(t_re) > 8.9884656743115785E+307) || (fabs(t_im) > 8.9884656743115785E+307))
			{
				b_t_re = t_re;
				t_re = log(rt_hypotd_snf(t_re / 2.0, t_im / 2.0)) + 0.69314718055994529;
				t_im = rt_atan2d_snf(t_im, b_t_re);
			}
			else
			{
				b_t_re = t_re;
				t_re = log(rt_hypotd_snf(t_re, t_im));
				t_im = rt_atan2d_snf(t_im, b_t_re);
			}

			zeta1->re = fnu * t_re - 0.0 * t_im;
			zeta1->im = fnu * t_im + 0.0 * t_re;
			zeta2->re = fnu * sr.re - 0.0 * sr.im;
			zeta2->im = fnu * sr.im + 0.0 * sr.re;

			if(sr.im == 0.0)
			{
				sr.re = rfn / sr.re;
				sr.im = 0.0;
			}
			else if(sr.re == 0.0)
			{
				if(rfn == 0.0)
				{
					sr.re = 0.0 / sr.im;
					sr.im = 0.0;
				}
				else
				{
					sr.re = 0.0;
					sr.im = -(rfn / sr.im);
				}
			}
			else
			{
				crfn_re = fabs(sr.re);
				b_t_re = fabs(sr.im);

				if(crfn_re > b_t_re)
				{
					ac = sr.im / sr.re;
					b_t_re = sr.re + ac * sr.im;
					sr.re = (rfn + ac * 0.0) / b_t_re;
					sr.im = (0.0 - ac * rfn) / b_t_re;
				}
				else if(b_t_re == crfn_re)
				{
					if(sr.re > 0.0)
						ac = 0.5;
					else
						ac = -0.5;

					if(sr.im > 0.0)
						b_t_re = 0.5;
					else
						b_t_re = -0.5;

					sr.re = (rfn * ac + 0.0 * b_t_re) / crfn_re;
					sr.im = (0.0 * ac - rfn * b_t_re) / crfn_re;
				}
				else
				{
					ac = sr.re / sr.im;
					b_t_re = sr.im + ac * sr.re;
					sr.re = ac * rfn / b_t_re;
					sr.im = (ac * 0.0 - rfn) / b_t_re;
				}
			}

			cwrk[15] = sr;
			b_sqrt(&cwrk[15]);
			phi->re = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) * cwrk[15].re;
			phi->im = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) * cwrk[15].im;

			if(ipmtr != 0)
				;
			else
			{
				if(s_im == 0.0)
				{
					cfn_re = 1.0 / s_re;
					t_im = 0.0;
				}
				else if(s_re == 0.0)
				{
					cfn_re = 0.0;
					t_im = -(1.0 / s_im);
				}
				else
				{
					crfn_re = fabs(s_re);
					b_t_re = fabs(s_im);

					if(crfn_re > b_t_re)
					{
						ac = s_im / s_re;
						b_t_re = s_re + ac * s_im;
						cfn_re = (1.0 + ac * 0.0) / b_t_re;
						t_im = (0.0 - ac) / b_t_re;
					}
					else if (b_t_re == crfn_re)
					{
						if(s_re > 0.0)
							ac = 0.5;
						else
							ac = -0.5;

						if(s_im > 0.0)
							b_t_re = 0.5;
						else
							b_t_re = -0.5;

						cfn_re = (ac + 0.0 * b_t_re) / crfn_re;
						t_im = (0.0 * ac - b_t_re) / crfn_re;
					}
					else
					{
						ac = s_re / s_im;
						b_t_re = s_im + ac * s_re;
						cfn_re = ac / b_t_re;
						t_im = (ac * 0.0 - 1.0) / b_t_re;
					}
				}

				cwrk[0].re = 1.0;
				cwrk[0].im = 0.0;
				crfn_re = 1.0;
				t_re = 0.0;
				ac = 1.0;
				l = 0;
				i = 1;

				do
				{
					exitg1 = 0;

					if(i < 15)
					{
						s_re = 0.0;
						s_im = 0.0;

						for(j = 0; j <= i; j++)
						{
							l++;
							b_t_re = s_re;
							s_re = (s_re * cfn_re - s_im * t_im) + C[l];
							s_im = b_t_re * t_im + s_im * cfn_re;
						}

						b_t_re = crfn_re;
						crfn_re = crfn_re * sr.re - t_re * sr.im;
						t_re = b_t_re * sr.im + t_re * sr.re;
						cwrk[i].re = crfn_re * s_re - t_re * s_im;
						cwrk[i].im = crfn_re * s_im + t_re * s_re;
						ac *= rfn;

						if((ac < tol) && (fabs(cwrk[i].re) + fabs(cwrk[i].im) < tol))
						{
							guard1 = true;
							exitg1 = 1;
						}
						else
							i++;
					}
					else
					{
						guard1 = true;
						exitg1 = 1;
					}
				}
				while(exitg1 == 0);
			}
		}
		else
		{
			zeta1->re = 1402.9773265065639 + fnu;
			zeta1->im = 0.0;
			zeta2->re = fnu;
			zeta2->im = 0.0;
			phi->re = 1.0;
			phi->im = 0.0;
		}
	}
	else
	{
		zeta1->re = 0.0;
		zeta1->im = 0.0;
		zeta2->re = 0.0;
		zeta2->im = 0.0;
		guard1 = true;
	}

	if(guard1)
	{
		if(ikflg == 2)
		{
			phi->re = 1.2533141373155003 * cwrk[15].re;
			phi->im = 1.2533141373155003 * cwrk[15].im;
		}
		else
		{
			phi->re = 0.3989422804014327 * cwrk[15].re;
			phi->im = 0.3989422804014327 * cwrk[15].im;
		}
	}
}

int b_cuoik(complex z, double fnu, int kode, int ikflg, int nin, complex y[2], double tol, double elim, double alim)
{
	int nuf;
	int n;
	complex zr;
	int iform;
	double gnn;
	double aarg;
	complex arg;
	complex an;
	complex phi;
	complex cz;
	complex zeta2;
	bool guard1 = false;
	complex ax;

	if(nin <= 2)
		n = nin;
	else
		n = 2;

	nuf = 0;

	if(z.re < 0.0)
	{
		zr.re = -z.re;
		zr.im = -z.im;
	}
	else
		zr = z;

	if(fabs(zr.im) > fabs(z.re) * 1.7321)
		iform = 2;
	else
		iform = 1;

	if(ikflg == 1)
	{
		if(fnu < 1.0)
			gnn = 1.0;
		else
			gnn = fnu;
	}
	else
	{
		gnn = (fnu + (double)n) - 1.0;

		if (gnn >= n)
			;
		else
			gnn = n;
	}

	aarg = 0.0;

	if(iform == 2)
	{
		if(zr.im <= 0.0)
		{
			an.re = -zr.im;
			an.im = -zr.re;
		}
		else
		{
			an.re = zr.im;
			an.im = -zr.re;
		}

		cunhj(an, gnn, tol, &phi, &arg, &cz, &zeta2);
		cz.re = zeta2.re - cz.re;
		cz.im = zeta2.im - cz.im;
		aarg = rt_hypotd_snf(arg.re, arg.im);
	}
	else
	{
		arg.re = 0.0;
		arg.im = 0.0;
		cunik(zr, gnn, ikflg, 1, tol, 0, &phi, &cz, &zeta2);
		cz.re = zeta2.re - cz.re;
		cz.im = zeta2.im - cz.im;
	}

	if(kode == 2)
	{
		cz.re -= zr.re;
		cz.im -= zr.im;
	}

	if(ikflg == 2)
	{
		cz.re = -cz.re;
		cz.im = -cz.im;
	}

	gnn = rt_hypotd_snf(phi.re, phi.im);

	if(cz.re >= alim)
	{
		gnn = cz.re + log(gnn);

		if(iform == 2)
			gnn = (gnn - 0.25 * log(aarg)) - 1.2655121234846454;

		if(gnn > elim)
			nuf = -1;
	}
	else
	{
		guard1 = false;

		if(cz.re >= -elim)
		{
			if(cz.re > -alim)
				;
			else
			{
				gnn = cz.re + log(gnn);

				if(iform == 2)
					gnn = (gnn - 0.25 * log(aarg)) - 1.2655121234846454;

				if(gnn > -elim)
				{
					if((phi.im == 0.0) && rtIsNaN(phi.re))
						;
					else if((fabs(phi.re) > 8.9884656743115785E+307) || (fabs(phi.im) > 8.9884656743115785E+307))
						phi.im = rt_atan2d_snf(phi.im, phi.re);
					else
						phi.im = rt_atan2d_snf(phi.im, phi.re);

					cz.im += phi.im;

					if(iform == 2)
					{
						phi = arg;

						if((arg.im == 0.0) && rtIsNaN(arg.re))
							;
						else if((fabs(arg.re) > 8.9884656743115785E+307) || (fabs(arg.im) > 8.9884656743115785E+307))
							phi.im = rt_atan2d_snf(arg.im, arg.re);
						else
							phi.im = rt_atan2d_snf(arg.im, arg.re);

						phi.im *= 0.25;
						phi.im = cz.im - phi.im;
						cz.im = phi.im;
					}

					gnn = exp(gnn) / tol;
					ax.re = gnn * cos(cz.im);
					ax.im = gnn * sin(cz.im);

					if(cuchk(ax, 2.2250738585072014E-305 / tol, tol) != 1)
						;
					else
						guard1 = true;
				}
				else
					guard1 = true;
			}
		}
		else
			guard1 = true;

		if(guard1)
		{
			for(iform = 1; iform <= n; iform++)
			{
				y[iform - 1].re = 0.0;
				y[iform - 1].im = 0.0;
			}

			nuf = n;
		}
	}

	return nuf;
}

int cuoik(complex z, double fnu, int kode, int ikflg, int nin, complex *y, double tol, double elim, double alim)
{
	int nuf;
	int n;
	complex zr;
	int iform;
	double gnn;
	double aarg;
	complex arg;
	complex an;
	complex phi;
	complex cz;
	complex zeta2;
	bool guard1 = false;
	complex ax;

	if(nin <= 1)
		n = nin;
	else
		n = 1;

	nuf = 0;

	if(z.re < 0.0)
	{
		zr.re = -z.re;
		zr.im = -z.im;
	}
	else
		zr = z;

	if(fabs(zr.im) > fabs(z.re) * 1.7321)
		iform = 2;
	else
		iform = 1;

	if(ikflg == 1)
	{
		if(fnu < 1.0)
			gnn = 1.0;
		else
			gnn = fnu;
	}
	else
	{
		gnn = (fnu + (double)n) - 1.0;

		if(gnn >= n)
			;
		else
			gnn = n;
	}

	aarg = 0.0;

	if(iform == 2)
	{
		if(zr.im <= 0.0)
		{
			an.re = -zr.im;
			an.im = -zr.re;
		}
		else
		{
			an.re = zr.im;
			an.im = -zr.re;
		}

		cunhj(an, gnn, tol, &phi, &arg, &cz, &zeta2);
		cz.re = zeta2.re - cz.re;
		cz.im = zeta2.im - cz.im;
		aarg = rt_hypotd_snf(arg.re, arg.im);
	}
	else
	{
		arg.re = 0.0;
		arg.im = 0.0;
		cunik(zr, gnn, ikflg, 1, tol, 0, &phi, &cz, &zeta2);
		cz.re = zeta2.re - cz.re;
		cz.im = zeta2.im - cz.im;
	}

	if(kode == 2)
	{
		cz.re -= zr.re;
		cz.im -= zr.im;
	}

	if(ikflg == 2)
	{
		cz.re = -cz.re;
		cz.im = -cz.im;
	}

	gnn = rt_hypotd_snf(phi.re, phi.im);

	if(cz.re >= alim)
	{
		gnn = cz.re + log(gnn);

		if(iform == 2)
			gnn = (gnn - 0.25 * log(aarg)) - 1.2655121234846454;

		if(gnn > elim)
			nuf = -1;
	}
	else
	{
		guard1 = false;

		if(cz.re >= -elim)
		{
			if(cz.re > -alim)
				;
			else
			{
				gnn = cz.re + log(gnn);

				if(iform == 2)
					gnn = (gnn - 0.25 * log(aarg)) - 1.2655121234846454;

				if(gnn > -elim)
				{
					if((phi.im == 0.0) && rtIsNaN(phi.re))
						;
					else if ((fabs(phi.re) > 8.9884656743115785E+307) || (fabs(phi.im) > 8.9884656743115785E+307))
						phi.im = rt_atan2d_snf(phi.im, phi.re);
					else
						phi.im = rt_atan2d_snf(phi.im, phi.re);

					cz.im += phi.im;

					if(iform == 2)
					{
						phi = arg;

						if((arg.im == 0.0) && rtIsNaN(arg.re))
							;
						else if((fabs(arg.re) > 8.9884656743115785E+307) || (fabs(arg.im) > 8.9884656743115785E+307))
							phi.im = rt_atan2d_snf(arg.im, arg.re);
						else
							phi.im = rt_atan2d_snf(arg.im, arg.re);

						phi.im *= 0.25;
						phi.im = cz.im - phi.im;
						cz.im = phi.im;
					}

					gnn = exp(gnn) / tol;
					ax.re = gnn * cos(cz.im);
					ax.im = gnn * sin(cz.im);

					if(cuchk(ax, 2.2250738585072014E-305 / tol, tol) != 1)
						;
					else
						guard1 = true;
				}
				else
					guard1 = true;
			}
		}
		else
			guard1 = true;

		if(guard1)
		{
			iform = 1;

			while(iform <= n)
			{
				y->re = 0.0;
				y->im = 0.0;
				iform = 2;
			}

			nuf = n;
		}
	}

	return nuf;
}

void b_exp(complex *x)
{
	double r;
	double x_im;

	if(rtIsInf(x->im) && rtIsInf(x->re) && (x->re < 0.0))
	{
		x->re = 0.0;
		x->im = 0.0;
	}
	else
	{
		r = exp(x->re / 2.0);
		x_im = x->im;
		x->re = r * (r * cos(x->im));
		x->im = r * (r * sin(x_im));
	}
}

void b_fix(double *x)
{
	if(*x < 0.0)
		*x = ceil(*x);
	else
		*x = floor(*x);
}

void gammaln(double *x)
{
	double t;
	double r;
	static const double table100[100] = { 0.0, 0.0, 0.69314718055994529,
		1.791759469228055, 3.1780538303479458, 4.7874917427820458,
		6.5792512120101012, 8.5251613610654147, 10.604602902745251,
		12.801827480081469, 15.104412573075516, 17.502307845873887,
		19.987214495661885, 22.552163853123425, 25.19122118273868, 27.89927138384089,
		30.671860106080672, 33.505073450136891, 36.395445208033053,
		39.339884187199495, 42.335616460753485, 45.380138898476908,
		48.471181351835227, 51.606675567764377, 54.784729398112319,
		58.003605222980518, 61.261701761002, 64.557538627006338, 67.88974313718154,
		71.257038967168015, 74.658236348830158, 78.0922235533153, 81.557959456115043,
		85.054467017581516, 88.580827542197682, 92.1361756036871, 95.7196945421432,
		99.330612454787428, 102.96819861451381, 106.63176026064346,
		110.32063971475739, 114.03421178146171, 117.77188139974507,
		121.53308151543864, 125.3172711493569, 129.12393363912722,
		132.95257503561632, 136.80272263732635, 140.67392364823425,
		144.5657439463449, 148.47776695177302, 152.40959258449735, 156.3608363030788,
		160.3311282166309, 164.32011226319517, 168.32744544842765,
		172.35279713916279, 176.39584840699735, 180.45629141754378,
		184.53382886144948, 188.6281734236716, 192.7390472878449, 196.86618167289,
		201.00931639928152, 205.1681994826412, 209.34258675253685,
		213.53224149456327, 217.73693411395422, 221.95644181913033,
		226.1905483237276, 230.43904356577696, 234.70172344281826,
		238.97838956183432, 243.26884900298271, 247.57291409618688,
		251.89040220972319, 256.22113555000954, 260.56494097186322,
		264.92164979855278, 269.29109765101981, 273.67312428569369,
		278.06757344036612, 282.4742926876304, 286.893133295427, 291.32395009427029,
		295.76660135076065, 300.22094864701415, 304.68685676566872,
		309.1641935801469, 313.65282994987905, 318.1526396202093, 322.66349912672615,
		327.1852877037752, 331.71788719692847, 336.26118197919845, 340.815058870799,
		345.37940706226686, 349.95411804077025, 354.53908551944079,
		359.1342053695754 };

	int i;
	static const double p1[8] = { 4.9452353592967269, 201.8112620856775,
		2290.8383738313464, 11319.672059033808, 28557.246356716354,
		38484.962284437934, 26377.487876241954, 7225.8139797002877 };

	static const double p2[8] = { 4.974607845568932, 542.4138599891071,
		15506.938649783649, 184793.29044456323, 1.0882047694688288E+6,
		3.33815296798703E+6, 5.1066616789273527E+6, 3.0741090548505397E+6 };

	static const double q1[8] = { 67.482125503037778, 1113.3323938571993,
		7738.7570569353984, 27639.870744033407, 54993.102062261576,
		61611.221800660023, 36351.2759150194, 8785.5363024310136 };

	static const double q2[8] = { 183.03283993705926, 7765.0493214450062,
		133190.38279660742, 1.1367058213219696E+6, 5.2679641174379466E+6,
		1.3467014543111017E+7, 1.7827365303532742E+7, 9.5330955918443538E+6 };

	static const double p4[8] = { 14745.0216605994, 2.4268133694867045E+6,
		1.2147555740450932E+8, 2.6634324496309772E+9, 2.9403789566345539E+10,
		1.7026657377653989E+11, 4.926125793377431E+11, 5.6062518562239514E+11 };

	static const double c[7] = { -0.001910444077728, 0.00084171387781295,
		-0.00059523799130430121, 0.0007936507935003503, -0.0027777777777776816,
		0.083333333333333329, 0.0057083835261 };

	static const double q4[8] = { 2690.5301758708993, 639388.56543000927,
		4.1355999302413881E+7, 1.120872109616148E+9, 1.4886137286788137E+10,
		1.0168035862724382E+11, 3.4174763455073773E+11, 4.4631581874197131E+11 };

	if(rtIsNaN(*x))
		;
	else if(*x < 0.0)
	{
		printf("NaN\n");
		*x = rtNaN;
	}
	else if(*x > 2.55E+305)
	{
		printf("Inf\n");
		*x = rtInf;
	}
	else if(*x <= 2.2204460492503131E-16)
	{
		printf("log\n");
		*x = -log(*x);
	}
	else if(*x <= 0.5)
	{
		printf("x <= 0.5\n");
		t = 1.0;
		r = 0.0;

		for(i = 0; i < 8; i++)
		{
			r = r * *x + p1[i];
			t = t * *x + q1[i];
		}

		*x = -log(*x) + *x * (-0.57721566490153287 + *x * (r / t));
	}
	else if(*x <= 0.6796875)
	{
		printf("x <= 0.6796875\n");
		t = 1.0;
		r = 0.0;

		for(i = 0; i < 8; i++)
		{
			r = r * ((*x - 0.5) - 0.5) + p2[i];
			t = t * ((*x - 0.5) - 0.5) + q2[i];
		}

		*x = -log(*x) + ((*x - 0.5) - 0.5) * (0.42278433509846713 + ((*x - 0.5) - 0.5) * (r / t));
	}
	else if((*x == floor(*x)) && (*x <= 100.0))
	{
		//printf("floor\n");
		*x = table100[(int)*x - 1];
	}
	else if(*x <= 1.5)
	{
		printf("x <= 1.5\n");
		t = 1.0;
		r = 0.0;

		for(i = 0; i < 8; i++)
		{
			r = r * ((*x - 0.5) - 0.5) + p1[i];
			t = t * ((*x - 0.5) - 0.5) + q1[i];
		}

		*x = ((*x - 0.5) - 0.5) * (-0.57721566490153287 + ((*x - 0.5) - 0.5) * (r / t));
	}
	else if (*x <= 4.0)
	{
		printf("x <= 4.0\n");
		t = 1.0;
		r = 0.0;

		for(i = 0; i < 8; i++)
		{
			r = r * (*x - 2.0) + p2[i];
			t = t * (*x - 2.0) + q2[i];
		}

		*x = (*x - 2.0) * (0.42278433509846713 + (*x - 2.0) * (r / t));
	}
	else if (*x <= 12.0)
	{
		printf("x <= 12.0\n");
		t = -1.0;
		r = 0.0;

		for(i = 0; i < 8; i++)
		{
			r = r * (*x - 4.0) + p4[i];
			t = t * (*x - 4.0) + q4[i];
		}

		*x = 1.791759469228055 + (*x - 4.0) * (r / t);
	}
	else
	{
		if(*x <= 2.25E+76)
		{
			printf("x <= 2.25E+76\n");
			r = 0.0057083835261;
			t = 1.0 / (*x * *x);

			for(i = 0; i < 6; i++)
				r = r * t + c[i];

			r /= *x;
		}
		else
			r = 0.0;

		t = log(*x);
		*x = ((r + 0.91893853320467278) - 0.5 * t) + *x * (t - 1.0);
	}
}

void b_log(complex *x)
{
	double x_im;
	double x_re;

	if((x->im == 0.0) && rtIsNaN(x->re))
	{
	}
	else if((fabs(x->re) > 8.9884656743115785E+307) || (fabs(x->im) > 8.9884656743115785E+307))
	{
		x_im = x->im;
		x_re = x->re;
		x->re = log(rt_hypotd_snf(x->re / 2.0, x->im / 2.0)) + 0.69314718055994529;
		x->im = rt_atan2d_snf(x_im, x_re);
	}
	else
	{
		x_im = x->im;
		x_re = x->re;
		x->re = log(rt_hypotd_snf(x->re, x->im));
		x->im = rt_atan2d_snf(x_im, x_re);
	}
}

double rtGetInf(void)
{
	size_t bitsPerReal = sizeof(double) * (NumBitsPerChar);
	double inf = 0.0;

	if(bitsPerReal == 32U)
		inf = rtGetInfF();
	else
	{
		unsigned short one = 1U;

		enum
		{
			LittleEndian,
			BigEndian
		}
		machByteOrder = (*((unsigned char *) &one) == 1U) ? LittleEndian : BigEndian;

		switch(machByteOrder)
		{
			case LittleEndian:
				{
					union
					{
						LittleEndianIEEEDouble bitVal;
						double fltVal;
					}
					tmpVal;

					tmpVal.bitVal.words.wordH = 0x7FF00000U;
					tmpVal.bitVal.words.wordL = 0x00000000U;
					inf = tmpVal.fltVal;
					break;
				}

			case BigEndian:
				{
					union
					{
						BigEndianIEEEDouble bitVal;
						double fltVal;
					}
					tmpVal;

					tmpVal.bitVal.words.wordH = 0x7FF00000U;
					tmpVal.bitVal.words.wordL = 0x00000000U;
					inf = tmpVal.fltVal;
					break;
				}
		}
	}

	return inf;
}

float rtGetInfF(void)
{
	IEEESingle infF;
	infF.wordL.wordLuint = 0x7F800000U;
	return infF.wordL.wordLreal;
}

double rtGetMinusInf(void)
{
	size_t bitsPerReal = sizeof(double) * (NumBitsPerChar);
	double minf = 0.0;

	if(bitsPerReal == 32U)
		minf = rtGetMinusInfF();
	else
	{
		unsigned short one = 1U;
		enum
		{
			LittleEndian,
			BigEndian
		}
		machByteOrder = (*((unsigned char *) &one) == 1U) ? LittleEndian : BigEndian;

		switch(machByteOrder)
		{
			case LittleEndian:
				{
					union
					{
						LittleEndianIEEEDouble bitVal;
						double fltVal;
					} tmpVal;

					tmpVal.bitVal.words.wordH = 0xFFF00000U;
					tmpVal.bitVal.words.wordL = 0x00000000U;
					minf = tmpVal.fltVal;
					break;
				}

			case BigEndian:
				{
					union
					{
						BigEndianIEEEDouble bitVal;
						double fltVal;
					}
					tmpVal;

					tmpVal.bitVal.words.wordH = 0xFFF00000U;
					tmpVal.bitVal.words.wordL = 0x00000000U;
					minf = tmpVal.fltVal;
					break;
				}
		}
	}

	return minf;
}

float rtGetMinusInfF(void)
{
	IEEESingle minfF;
	minfF.wordL.wordLuint = 0xFF800000U;
	return minfF.wordL.wordLreal;
}

double rtGetNaN(void)
{
	size_t bitsPerReal = sizeof(double) * (NumBitsPerChar);
	double nan = 0.0;

	if(bitsPerReal == 32U)
		nan = rtGetNaNF();
	else
	{
		unsigned short one = 1U;
		enum
		{
			LittleEndian,
			BigEndian
		}
		machByteOrder = (*((unsigned char *) &one) == 1U) ? LittleEndian : BigEndian;

		switch(machByteOrder)
		{
			case LittleEndian:
				{
					union
					{
						LittleEndianIEEEDouble bitVal;
						double fltVal;
					}
					tmpVal;

					tmpVal.bitVal.words.wordH = 0xFFF80000U;
					tmpVal.bitVal.words.wordL = 0x00000000U;
					nan = tmpVal.fltVal;
					break;
				}

			case BigEndian:
				{
					union
					{
						BigEndianIEEEDouble bitVal;
						double fltVal;
					}
					tmpVal;

					tmpVal.bitVal.words.wordH = 0x7FFFFFFFU;
					tmpVal.bitVal.words.wordL = 0xFFFFFFFFU;
					nan = tmpVal.fltVal;
					break;
				}
		}
	}

	return nan;
}

float rtGetNaNF(void)
{
	IEEESingle nanF = { { 0 } };

	unsigned short one = 1U;
	enum
	{
		LittleEndian,
		BigEndian
	}
	machByteOrder = (*((unsigned char *) &one) == 1U) ? LittleEndian : BigEndian;

	switch(machByteOrder)
	{
		case LittleEndian:
			{
				nanF.wordL.wordLuint = 0xFFC00000U;
				break;
			}

		case BigEndian:
			{
				nanF.wordL.wordLuint = 0x7FFFFFFFU;
				break;
			}
	}

	return nanF.wordL.wordLreal;
}

void rt_InitInfAndNaN(size_t realSize)
{
	(void) (realSize);
	rtNaN = rtGetNaN();
	rtNaNF = rtGetNaNF();
	rtInf = rtGetInf();
	rtInfF = rtGetInfF();
	rtMinusInf = rtGetMinusInf();
	rtMinusInfF = rtGetMinusInfF();
}

bool rtIsInf(double value)
{
	return ((value==rtInf || value==rtMinusInf) ? 1U : 0U);
}

bool rtIsInfF(float value)
{
	return (((value)==rtInfF || (value)==rtMinusInfF) ? 1U : 0U);
}

bool rtIsNaN(double value)
{
	return (value!=value)? 1U:0U;
}

bool rtIsNaNF(float value)
{
	return (value!=value)? 1U:0U;
}

void b_sqrt(complex *x)
{
	double absxi;
	double absxr;

	if(x->im == 0.0)
	{
		printf("x.re == 0.0\n");

		if(x->re < 0.0)
		{
			printf("x.re < 0.0\n");
			absxi = 0.0;
			absxr = sqrt(fabs(x->re));
		}
		else
		{
			printf("else x.re < 0.0\n");
			absxi = sqrt(x->re);
			absxr = 0.0;
		}
	}
	else if(x->re == 0.0)
	{
		printf("x.re == 0.0\n");

		if(x->im < 0.0)
		{
			printf("x.im < 0.0\n");
			absxi = sqrt(-x->im / 2.0);
			absxr = -absxi;
		}
		else
		{
			printf("else x.im < 0.0\n");
			absxi = sqrt(x->im / 2.0);
			absxr = absxi;
		}
	}
	else if(rtIsNaN(x->re) || rtIsNaN(x->im))
	{
		absxi = rtNaN;
		absxr = rtNaN;
	}
	else if(rtIsInf(x->im))
	{
		absxi = rtInf;
		absxr = x->im;
	}
	else if(rtIsInf(x->re))
	{
		if(x->re < 0.0)
		{
			absxi = 0.0;
			absxr = rtInf;
		}
		else
		{
			absxi = rtInf;
			absxr = 0.0;
		}
	}
	else
	{
		printf("else\n");
		absxr = fabs(x->re);
		absxi = fabs(x->im);

		if((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307))
		{
			absxr *= 0.5;
			absxi *= 0.5;
			absxi = rt_hypotd_snf(absxr, absxi);

			if(absxi > absxr)
				absxi = sqrt(absxi) * sqrt(1.0 + absxr / absxi);
			else
				absxi = sqrt(absxi) * 1.4142135623730951;
		}
		else
			absxi = sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);

		if(x->re > 0.0)
			absxr = 0.5 * (x->im / absxi);
		else
		{
			if(x->im < 0.0)
				absxr = -absxi;
			else
				absxr = absxi;

			absxi = 0.5 * (x->im / absxr);
		}
	}

	x->re = absxi;
	x->im = absxr;
}

double rt_atan2d_snf(double u0, double u1)
{
	double y;
	int b_u0;
	int b_u1;

	if(rtIsNaN(u0) || rtIsNaN(u1))
		y = rtNaN;
	else if(rtIsInf(u0) && rtIsInf(u1))
	{
		if (u0 > 0.0)
			b_u0 = 1;
		else
			b_u0 = -1;

		if (u1 > 0.0)
			b_u1 = 1;
		else
			b_u1 = -1;

		y = atan2(b_u0, b_u1);
	}
	else if(u1 == 0.0)
	{
		if(u0 > 0.0)
			y = RT_PI / 2.0;
		else if(u0 < 0.0)
			y = -(RT_PI / 2.0);
		else
			y = 0.0;
	}
	else
		y = atan2(u0, u1);

	return y;
}

double rt_hypotd_snf(double u0, double u1)
{
	double y;
	double a;
	double b;
	a = fabs(u0);
	b = fabs(u1);

	if(a < b)
	{
		a /= b;
		y = b * sqrt(a * a + 1.0);
	}
	else if(a > b)
	{
		b /= a;
		y = a * sqrt(b * b + 1.0);
	}
	else if(rtIsNaN(b))
		y = b;
	else
		y = a * 1.4142135623730951;

	return y;
}

int cacai(complex z, double fnu, int kode, int mr, complex *y, double rl, double tol, double elim, double alim)
{
	int nz;
	complex zn;
	double az;
	bool guard1 = false;
	int nw;
	double crsc_re;
	bool iflag;
	double hz_re;
	double hz_im;
	double cz_re;
	complex cy[2];
	double cz_im;
	double acz;
	double ck_re;
	double aa;
	double ak1_re;
	double ak1_im;
	double ascle;
	complex s2;
	double coef_im;
	double b_atol;
	int exitg1;
	double s;
	double rs;
	double b_ak1_re;
	nz = 0;
	zn.re = -z.re;
	zn.im = -z.im;
	az = rt_hypotd_snf(z.re, z.im);
	guard1 = false;

	if((!(az <= 2.0)) && (az * az * 0.25 > ((fnu + 1.0) - 1.0) + 1.0))
	{
		if(az < rl)
			nw = cmlri(zn, fnu, kode, 1, y, tol);
		else
			nw = casyi(zn, fnu, kode, 1, y, rl, tol, elim);

		if(nw < 0)
		{
			if(nw == -2)
				nz = -2;
			else
				nz = -1;
		}
		else
			guard1 = true;
	}
	else
	{
		az = rt_hypotd_snf(-z.re, -z.im);

		if(az == 0.0)
		{
			if(fnu == 0.0)
			{
				y->re = 1.0;
				y->im = 0.0;
			}
			else
			{
				y->re = 0.0;
				y->im = 0.0;
			}

			guard1 = true;
		}
		else
		{
			crsc_re = 1.0;
			iflag = false;

			if(az < 2.2250738585072014E-305)
			{
				if(fnu == 0.0)
				{
					y->re = 1.0;
					y->im = 0.0;
				}
				else
				{
					y->re = 0.0;
					y->im = 0.0;
				}

				guard1 = true;
			}
			else
			{
				hz_re = 0.5 * -z.re;
				hz_im = 0.5 * -z.im;

				if(az > 4.7170688552396617E-153)
				{
					cz_re = hz_re * hz_re - hz_im * hz_im;
					cz_im = hz_re * hz_im + hz_im * hz_re;
					acz = rt_hypotd_snf(cz_re, cz_im);
				}
				else
				{
					cz_re = 0.0;
					cz_im = 0.0;
					acz = 0.0;
				}

				ck_re = hz_re;

				if((hz_im == 0.0) && rtIsNaN(hz_re))
					;
				else if((fabs(hz_re) > 8.9884656743115785E+307) || (fabs(hz_im) > 8.9884656743115785E+307))
				{
					ck_re = log(rt_hypotd_snf(hz_re / 2.0, hz_im / 2.0)) + 0.69314718055994529;
					hz_im = rt_atan2d_snf(hz_im, hz_re);
				}
				else
				{
					ck_re = log(rt_hypotd_snf(hz_re, hz_im));
					hz_im = rt_atan2d_snf(hz_im, hz_re);
				}

				az = ((fnu + 1.0) - 1.0) + 1.0;
				gammaln(&az);
				ak1_re = ck_re * ((fnu + 1.0) - 1.0) - az;
				ak1_im = hz_im * ((fnu + 1.0) - 1.0);

				if(kode == 2)
					ak1_re -= -z.re;

				if(ak1_re > -elim)
				{
					ascle = 0.0;

					if(ak1_re <= -alim)
					{
						iflag = true;
						crsc_re = tol;
						ascle = 2.2250738585072014E-305 / tol;
					}

					aa = exp(ak1_re);

					if(iflag)
						aa /= tol;

					hz_re = aa * cos(ak1_im);
					coef_im = aa * sin(ak1_im);
					b_atol = tol * acz / (((fnu + 1.0) - 1.0) + 1.0);
					nw = 1;

					do
					{
						exitg1 = 0;

						if(nw <= 1)
						{
							ck_re = 1.0;
							hz_im = 0.0;

							if(!(acz < tol * (fnu + 1.0)))
							{
								ak1_re = 1.0;
								ak1_im = 0.0;
								az = (fnu + 1.0) + 2.0;
								s = fnu + 1.0;
								aa = 2.0;

								do
								{
									rs = 1.0 / s;
									b_ak1_re = ak1_re;
									ak1_re = rs * (ak1_re * cz_re - ak1_im * cz_im);
									ak1_im = rs * (b_ak1_re * cz_im + ak1_im * cz_re);
									ck_re += ak1_re;
									hz_im += ak1_im;
									s += az;
									az += 2.0;
									aa = aa * acz * rs;
								}
								while(!!(aa > b_atol));
							}

							s2.re = ck_re * hz_re - hz_im * coef_im;
							s2.im = ck_re * coef_im + hz_im * hz_re;

							if(iflag && (cuchk(s2, ascle, tol) != 0))
							{
								guard1 = true;
								exitg1 = 1;
							}
							else
							{
								y->re = s2.re * crsc_re - s2.im * 0.0;
								y->im = s2.re * 0.0 + s2.im * crsc_re;
								nw = 2;
							}
						}
						else
						{
							guard1 = true;
							exitg1 = 1;
						}
					}
					while(exitg1 == 0);
				}
				else
				{
					y->re = 0.0;
					y->im = 0.0;
					guard1 = true;
				}
			}
		}
	}

	if(guard1)
	{
		for(nw = 0; nw < 2; nw++)
		{
			cy[nw].re = 0.0;
			cy[nw].im = 0.0;
		}

		nw = b_cbknu(zn, fnu, kode, 1, cy, tol, elim, alim);

		if(nw != 0)
		{
			if(nw == -2)
				nz = -2;
			else
				nz = -1;
		}
		else
		{
			if(mr < 0)
				aa = 3.1415926535897931;
			else
				aa = -3.1415926535897931;

			ck_re = 0.0;
			hz_im = aa;

			if(kode != 1)
			{
				az = cos(-(-z.im));
				hz_re = sin(-(-z.im));
				ck_re = 0.0 * az - aa * hz_re;
				hz_im = 0.0 * hz_re + aa * az;
			}

			nw = (int)fnu;
			az = (fnu - (double)nw) * aa;
			s2.re = cos(az);
			s2.im = sin(az);

			if((nw & 1) != 0)
			{
				s2.re = -s2.re;
				s2.im = -s2.im;
			}

			ak1_re = cy[0].re;
			ak1_im = cy[0].im;

			if(kode != 1)
			{
				ascle = 2.2250738585072014E-305 / tol;
				ak1_re = cy[0].re;
				ak1_im = cy[0].im;
				az = rt_hypotd_snf(cy[0].re, cy[0].im);

				if(az > 0.0)
				{
					if((-(-z.re) - (-z.re)) + log(az) < -alim)
					{
						ak1_re = 0.0;
						ak1_im = 0.0;
						az = 0.0;
					}
					else
					{
						if((cy[0].im == 0.0) && rtIsNaN(cy[0].re))
							;
						else if((fabs(cy[0].re) > 8.9884656743115785E+307) || (fabs(cy[0].im) > 8.9884656743115785E+307))
						{
							ak1_re = log(rt_hypotd_snf(cy[0].re / 2.0, cy[0].im / 2.0)) + 0.69314718055994529;
							ak1_im = rt_atan2d_snf(cy[0].im, cy[0].re);
						}
						else
						{
							ak1_re = log(rt_hypotd_snf(cy[0].re, cy[0].im));
							ak1_im = rt_atan2d_snf(cy[0].im, cy[0].re);
						}

						ak1_re -= -z.re;
						ak1_im -= -z.im;
						ak1_re -= -z.re;
						ak1_im -= -z.im;

						if(rtIsInf(ak1_im) && rtIsInf(ak1_re) && (ak1_re < 0.0))
						{
							ak1_re = 0.0;
							ak1_im = 0.0;
						}
						else
						{
							az = exp(ak1_re / 2.0);
							ak1_re = az * (az * cos(ak1_im));
							ak1_im = az * (az * sin(ak1_im));
						}

						az = rt_hypotd_snf(ak1_re, ak1_im);
					}
				}

				if((az > ascle) || (rt_hypotd_snf(y->re, y->im) > ascle))
					nz = 0;
				else
				{
					ak1_re = 0.0;
					ak1_im = 0.0;
					y->re = 0.0;
					y->im = 0.0;
					nz = 1;
				}
			}

			az = ck_re * y->im + hz_im * y->re;
			y->re = (ck_re * y->re - hz_im * y->im) + (s2.re * ak1_re - s2.im * ak1_im);
			y->im = az + (s2.re * ak1_im + s2.im * ak1_re);
		}
	}

	return nz;
}

complex cairy(const complex z, int id, int kode)
{
	complex ai;
	double az;
	double s1_re;
	double az3;
	double s1_im;
	double r;
	complex s2;
	complex trm2;
	double aa;
	int iflag;
	complex trm1;
	double ak;
	double atrm;
	bool guard1 = false;
	double z3_re;
	bool guard2 = false;
	double z3_im;
	bool guard3 = false;
	double d2;
	double ad;
	int nn;
	double bk;
	bool exitg1;
	int b_z;
	double b_z3_re;
	double b_z3_im;
	ai.re = 0.0;
	ai.im = 0.0;
	az = rt_hypotd_snf(z.re, z.im);

	if(az > 1.0)
	{
		az3 = (1.0 + (double)id) / 3.0;
		r = log(az);
		trm2 = z;
		b_sqrt(&trm2);
		s2.re = 0.66666666666666663 * (z.re * trm2.re - z.im * trm2.im);
		s2.im = 0.66666666666666663 * (z.re * trm2.im + z.im * trm2.re);
		iflag = 0;
		aa = 1.0;
		ak = s2.im;

		if(z.re < 0.0)
			s2.re = -fabs(s2.re);

		if((z.im != 0.0) || (z.re > 0.0))
			;
		else
		{
			s2.re = 0.0;
			s2.im = ak;
		}

		guard1 = false;
		guard2 = false;
		guard3 = false;

		if((s2.re >= 0.0) && (z.re > 0.0))
		{
			if((kode != 2) && (s2.re >= 664.87164553371019))
			{
				iflag = 2;
				aa = 4.503599627370496E+15;

				if(-s2.re - 0.25 * r < -700.92179369444591)
					;
				else
					guard3 = true;
			}
			else
				guard3 = true;
		}
		else if((kode != 2) && (s2.re <= -664.87164553371019))
		{
			iflag = 1;
			aa = 2.2204460492503131E-16;

			if(-s2.re + 0.25 * r > 700.92179369444591)
				;
			else
				guard2 = true;
		}
		else
			guard2 = true;

		if(guard3)
		{
			cbknu(s2, az3, kode, 664.87164553371019, &trm1, &nn);
			guard1 = true;
		}

		if(guard2)
		{
			trm1.re = 0.0;
			trm1.im = 0.0;

			if(z.im < 0.0)
				b_z = -1;
			else
				b_z = 1;

			nn = cacai(s2, az3, kode, b_z, &trm1, 21.784271729432426,
					2.2204460492503131E-16, 700.92179369444591, 664.87164553371019);
			if(nn < 0)
				;
			else
				guard1 = true;
		}

		if(guard1)
		{
			s1_re = 0.18377629847393068 * trm1.re;
			s1_im = 0.18377629847393068 * trm1.im;

			if(iflag != 0)
			{
				s1_re *= aa;
				s1_im *= aa;

				if(id == 1)
				{
					s1_re = -s1_re;
					s1_im = -s1_im;
					r = s1_re;
					s1_re = s1_re * z.re - s1_im * z.im;
					s1_im = r * z.im + s1_im * z.re;
				}
				else
				{
					r = s1_re;
					s1_re = s1_re * trm2.re - s1_im * trm2.im;
					s1_im = r * trm2.im + s1_im * trm2.re;
				}

				ai.re = 1.0 / aa * s1_re;
				ai.im = 1.0 / aa * s1_im;
			}
			else if(id == 1)
			{
				ai.re = -z.re * s1_re - -z.im * s1_im;
				ai.im = -z.re * s1_im + -z.im * s1_re;
			}
			else
			{
				ai.re = trm2.re * s1_re - trm2.im * s1_im;
				ai.im = trm2.re * s1_im + trm2.im * s1_re;
			}
		}
	}
	else
	{
		s1_re = 1.0;
		s1_im = 0.0;
		s2.re = 1.0;
		s2.im = 0.0;

		if(az < 2.2204460492503131E-16)
		{
			s1_re = 0.0;

			if(id == 1)
			{
				if(az > 4.7170688552396617E-153)
				{
					s1_re = 0.5 * (z.re * z.re - z.im * z.im);
					s1_im = 0.5 * (z.re * z.im + z.im * z.re);
				}

				s1_re *= 0.35502805388781722;
				s1_im *= 0.35502805388781722;
				ai.re = -0.25881940379280682 + s1_re;
				ai.im = s1_im;
			}
			else
			{
				if(az > 2.2250738585072014E-305)
				{
					s1_re = 0.25881940379280682 * z.re;
					s1_im = 0.25881940379280682 * z.im;
				}

				ai.re = 0.35502805388781722 - s1_re;
				ai.im = 0.0 - s1_im;
			}
		}
		else
		{
			aa = az * az;

			if(aa >= 2.2204460492503131E-16 / az)
			{
				trm1.re = 1.0;
				trm1.im = 0.0;
				trm2.re = 1.0;
				trm2.im = 0.0;
				atrm = 1.0;
				r = z.re * z.re - z.im * z.im;
				az3 = z.re * z.im + z.im * z.re;
				z3_re = r * z.re - az3 * z.im;
				z3_im = r * z.im + az3 * z.re;
				az3 = az * aa;
				aa = (2.0 + (double)id) * ((3.0 + (double)id) + (double)id);
				d2 = ((3.0 - (double)id) - (double)id) * (4.0 - (double)id);

				if(aa <= d2)
					ad = aa;
				else
					ad = d2;

				ak = 24.0 + 9.0 * (double)id;
				bk = 30.0 - 9.0 * (double)id;
				nn = 1;
				exitg1 = false;

				while((!exitg1) && (nn < 26))
				{
					if(z3_im == 0.0)
					{
						b_z3_re = z3_re / aa;
						b_z3_im = 0.0;
					}
					else if(z3_re == 0.0)
					{
						b_z3_re = 0.0;
						b_z3_im = z3_im / aa;
					}
					else
					{
						b_z3_re = z3_re / aa;
						b_z3_im = z3_im / aa;
					}

					r = trm1.re;
					trm1.re = trm1.re * b_z3_re - trm1.im * b_z3_im;
					trm1.im = r * b_z3_im + trm1.im * b_z3_re;
					s1_re += trm1.re;
					s1_im += trm1.im;

					if(z3_im == 0.0)
					{
						b_z3_re = z3_re / d2;
						b_z3_im = 0.0;
					}
					else if(z3_re == 0.0)
					{
						b_z3_re = 0.0;
						b_z3_im = z3_im / d2;
					}
					else
					{
						b_z3_re = z3_re / d2;
						b_z3_im = z3_im / d2;
					}

					r = trm2.re;
					trm2.re = trm2.re * b_z3_re - trm2.im * b_z3_im;
					trm2.im = r * b_z3_im + trm2.im * b_z3_re;
					s2.re += trm2.re;
					s2.im += trm2.im;
					atrm = atrm * az3 / ad;
					aa += ak;
					d2 += bk;

					if((aa <= d2) || rtIsNaN(d2))
						ad = aa;
					else
						ad = d2;

					if(atrm < 2.2204460492503131E-16 * ad)
						exitg1 = true;
					else
					{
						ak += 18.0;
						bk += 18.0;
						nn++;
					}
				}
			}

			if(id == 1)
			{
				ai.re = -0.25881940379280682 * s2.re;
				ai.im = -0.25881940379280682 * s2.im;

				if(az > 2.2204460492503131E-16)
				{
					r = z.re * z.re - z.im * z.im;
					az3 = z.re * z.im + z.im * z.re;
					ai.re += (r * s1_re - az3 * s1_im) * 0.17751402694390861;
					ai.im += (r * s1_im + az3 * s1_re) * 0.17751402694390861;
				}

				if(kode == 1)
					;
				else
				{
					s2 = z;
					b_sqrt(&s2);
					r = s2.re * z.im + s2.im * z.re;
					s2.re = 0.66666666666666663 * (s2.re * z.re - s2.im * z.im);
					s2.im = 0.66666666666666663 * r;

					if(rtIsInf(s2.im) && rtIsInf(s2.re) && (s2.re < 0.0))
					{
						s2.re = 0.0;
						s2.im = 0.0;
					}
					else
					{
						r = exp(s2.re / 2.0);
						s2.re = r * (r * cos(s2.im));
						s2.im = r * (r * sin(s2.im));
					}

					r = ai.re;
					ai.re = ai.re * s2.re - ai.im * s2.im;
					ai.im = r * s2.im + ai.im * s2.re;
				}
			}
			else
			{
				ai.re = s1_re * 0.35502805388781722 - (z.re * s2.re - z.im * s2.im) * 0.25881940379280682;
				ai.im = s1_im * 0.35502805388781722 - (z.re * s2.im + z.im * s2.re) * 0.25881940379280682;

				if(kode == 1)
					;
				else
				{
					s2 = z;
					b_sqrt(&s2);
					r = s2.re * z.im + s2.im * z.re;
					s2.re = 0.66666666666666663 * (s2.re * z.re - s2.im * z.im);
					s2.im = 0.66666666666666663 * r;

					if(rtIsInf(s2.im) && rtIsInf(s2.re) && (s2.re < 0.0))
					{
						s2.re = 0.0;
						s2.im = 0.0;
					}
					else
					{
						r = exp(s2.re / 2.0);
						s2.re = r * (r * cos(s2.im));
						s2.im = r * (r * sin(s2.im));
					}

					r = ai.re;
					ai.re = ai.re * s2.re - ai.im * s2.im;
					ai.im = r * s2.im + ai.im * s2.re;
				}
			}
		}
	}

	return ai;
}
