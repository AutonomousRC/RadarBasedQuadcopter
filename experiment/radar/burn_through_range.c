#include "radar_basic.h"
#include "math_tech.h"

void burn_through_range(float pt, float g, float sigma, float freq, float tau, float loss, float pj,
                        float bj, float gj, float lossj, float sdjpn, float erp, float *range)
{
	float k, sjn, lambda, gje, ge, Ar, num32, demo3, demo4, val1, val2, val3, val4;

	k = 1.38 * pow(10, -23);
	sjn = pow(10, sdjpn / 10);
	lambda = get_wavelength(freq);
	gje = pow(10, gj / 10);
	ge = pow(10, g / 10);
	Ar = pow(lambda, 2) * ge / (4 * M_PI);
	num32 = erp * Ar;
	demo3 = 8 * M_PI * bj * k * 290;
	demo4 = pow(4 * M_PI, 2) * k * 290 * sjn;
	val1 = pow((num32 / demo3), 2);
	val2 = (pt * tau * ge * sigma * Ar) / (pow((4 * M_PI), 2) * loss * sjn * k * 290);
	val3 = sqrt(val1 + val2);
	val4 = (erp * Ar) / demo3;

	*range = sqrt(val3 - val4) / 1000;
}
