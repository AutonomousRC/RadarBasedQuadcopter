#include "init_require_data.h"
#include "bluestein_fft.h"
#include "complex_math.h"
#include "ogl_helper.h"
#include "math_tech.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void arr_alloc(double **data, int num)
{
	*data = (double *)malloc(sizeof(double) * num);
	memset(*data, 0x0, sizeof(double) * num);
}

void complex_arr_alloc(comp **data, int num)
{
	*data = (comp *)malloc(sizeof(comp) * num);
	memset(*data, 0x0, sizeof(comp) * num);
}

// first sample_num = 800
void init_require_data(double *sig, int sample_num)
{
	int i, n = sample_num;

	double tau = 10 * pow(10, -6);
	double b = 40 * pow(10, 6);
	double time_b_prod = b * tau;
	//double t[n] = {0};
	double *t = NULL;
	//double sq_t[n] = {0};
	double *sq_t = NULL;
	//double series[n] = {0};
	double *series = NULL;

	//comp data[n] = {0};
	comp *data = NULL;
	//comp fft_res[n] = {0};
	comp *fft_res = NULL;

	arr_alloc(&t, n);
	arr_alloc(&sq_t, n);
	arr_alloc(&series, n);

	complex_arr_alloc(&data, n);
	complex_arr_alloc(&fft_res, n);

	linear_slice(-tau / 2.0, tau / 2.0, n, t);
	//p_arr(t, n);
	pow_vec(t, sq_t, n, 2.0);

	for(i = 0; i < n; i++)
		series[i] = sq_t[i] * M_PI * (b / tau);

	//p_arr(series, n);

	euler_formula(series, data, n);

	//print_complex(data, n);

	bluestein_fft(data, fft_res, n);

	/* This */
	//print_complex(fft_res, n);

	complex_abs(fft_res, sig, n);

	//p_arr(sig, n);

	fft_shift(sig);

	p_arr(sig, n);
}
