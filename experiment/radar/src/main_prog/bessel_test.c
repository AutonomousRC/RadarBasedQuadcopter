#include "zeroth_modify_bessel.h"
#include "complex_math.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
	double data[33] = {0};
	int i, len = 32;

	double tmp1[33] = {0};
	double tmp2[33] = {0};
	double tmp3[33] = {0};
	double tmp4[33] = {0};
	double tmp5[33] = {0};
	complex tmp6[33] = {0};
	complex z = {5.0, 0.0};

	for(i = 0; i < len; i++)
	{
		tmp1[i] = i - (len - 1) / 2.0;
		tmp2[i] = (len - 1) / 2.0;
		tmp3[i] = tmp1[i] / tmp2[i];
		tmp4[i] = pow(tmp3[i], 2.0);
		tmp5[i] = 5 * sqrt(1 - tmp4[i]);
		//printf("tmp1[%d] = %lf, tmp2[%d] = %lf\n", i, tmp1[i], i, tmp2[i]);
		//printf("tmp3[%d] = %lf\n", i, tmp3[i]);
		//printf("tmp4[%d] = %lf\n", i, tmp4[i]);
		//printf("tmp5[%d] = %lf\n", i, tmp5[i]);
		//print_complex(tmp6, 32);
	}

	zeroth_modify_bessel(tmp5, tmp6, 32);

#if 0
	for(i = 0; i < len; i++)
		print_complex(tmp6, 32);
#endif
		

	return 0;
}
