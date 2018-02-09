#include "init_require_data.h"
#include "complex_math.h"
#include "math_tech.h"
#include <math.h>

void init_require_data(double *sig)
{
	int i, n = 800;

	double tau = 10 * pow(10, -6);
	double b = 40 * pow(10, 6);
	double time_b_prod = b * tau;
	double t[800] = {0};
	double sq_t[800] = {0};
	double series[800] = {0};

	complex data[800] = {0};
	complex fft_res[800] = {0};
	
	linear_slice(-tau / 2.0, tau / 2.0, n, t);
	pow_vec(t, sq_t, 800, 2.0);

	for(i = 0; i < 800; i++)
		series[i] = sq_t[i] * M_PI * (b / tau);

	euler_formula(series, data, 800);

	//print_complex(data, 800);

	fft(data, fft_res);
}
