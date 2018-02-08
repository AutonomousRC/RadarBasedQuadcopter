#include "radar_basic.h"
#include "math_tech.h"

void sdjpn(float pt, float g, float sigma, float freq, float tau, float loss,
	float R, float pj, float bj, float gj, float lossj, float *sjn)
{
	float k, range, lambda, gje, ge, ERP, ERP_db, Ar, num1, demo1, demo2, num2, val11, val21;

	k = 1.38 * pow(10, -23);
	range = R * 1000;
	lambda = get_wavelength(freq);
	gje = pow(10, gj / 10);
	ge = pow(10, g / 10);
	ERP = pj * gje / lossj;
	ERP_db = 10 * log10(ERP);
	Ar = pow(lambda, 2) * ge / (4 * M_PI);
	num1 = pt * tau * ge * sigma * Ar;
	demo1 = pow(4, 2) * pow(M_PI, 2) * loss * pow(range, 4);
	demo2 = 4 * M_PI * bj * pow(range, 2);
	num2 = ERP * Ar;
	val11 = num1 / demo1;
	val21 = num2 / demo2;
	*sjn = 10 * log10(val11 / (val21 + k * 290));
}
