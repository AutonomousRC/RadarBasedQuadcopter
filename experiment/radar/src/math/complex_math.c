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

void print_complex(complex *data, int data_num)
{
	int i;

	for(i = 0; i < data_num; i++)
	{
		printf("real = %lf\n", data[i].re);
		printf("imag = %lf\n", data[i].im);
	}
}
