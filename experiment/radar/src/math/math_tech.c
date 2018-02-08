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
