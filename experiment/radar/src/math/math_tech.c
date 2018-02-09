#include <math.h>
#include "math_tech.h"
#include "radar_basic.h"

float angle2radian(float ang)
{
        return ang * M_PI / 180.0;
}

float get_wavelength(float freq)
{
        return LIGHT_SPEED / freq;
}

void slice_section(double s, double e, double itv, double *arr)
{
        int i;
        double tmp = e / itv;

        for(i = 0; i < (int)(tmp + 2); i++)
        {
                arr[i] = itv * i;
#if 0
                printf("arr[%d] = %.9lf\n", i, arr[i]);
#endif
        }
}

void linear_slice(double start, double end, int num, double *arr)
{
        int i;
        double interval = (end - start) / (num - 1);
	//double interval = (end - start) / num;

        for(i = 0; i < num; i++)
                arr[i] = start + interval * i;
}

void pow_vec(double *x, double *y, int data_num, double series)
{
	int i;

	for(i = 0; i < data_num; i++)
		y[i] = pow(x[i], series);
}
