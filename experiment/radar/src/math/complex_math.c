#include <math.h>
#include <stdio.h>
#include "complex_math.h"

void euler_formula(double *series, complex *data, int data_num)
{
	int i;

	for(i = 0; i < data_num; i++)
	{
		data[i].re = cos(series[i]);
		data[i].im = sin(series[i]);
	}
}

void complex_abs(complex *x, double *y, int data_num)
{
	int i;
	double temp_re, temp_im;

	//for(i = 0; i < 800; i++)
	for(i = 0; i < data_num; i++)
	{
		temp_re = fabs(x[i].re);
		temp_im = fabs(x[i].im);

		y[i] = sqrt(pow(temp_re, 2) + pow(temp_im, 2));
	}
}

void print_complex(complex *data, int data_num)
{
	int i;

	for(i = 0; i < data_num; i++)
	{
		printf("real = %lf\n", data[i].re);
		printf("imag = %lf\n", data[i].im);
	}
}
