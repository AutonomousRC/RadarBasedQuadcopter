#include "radar_basic.h"
#include "math_tech.h"

void stand_off_jammer(float pt, float g, float sigma, float b, float freq, float loss, float range,
                float pj, float bj, float gj, float lossj, float gprime, float rangej, float *br_range)
{
	float lambda_sqdb, sigma_db, range_db, rangej_db, pt_db, b_db, bj_db, pj_db, factor;

	lambda_sqdb = 10 * log10(pow(get_wavelength(freq), 2));

	if(loss == 0.0)
		loss = 0.000001;

	if(lossj == 0.0)
		lossj = 0.000001;

	sigma_db = 10 * log10(sigma);
	range_db = 10 * log10(range * 1000.0);
	rangej_db = 10 * log10(rangej * 1000.0);
	pt_db = 10 * log10(pt);
	b_db = 10 * log10(b);
	bj_db = 10 * log10(bj);
	pj_db = 10 * log10(pj);
	factor = 10 * log10(4.0 * M_PI);

	*br_range = pow((pt * pow(10, 2 * g / 10) * sigma * bj * pow(10, lossj / 10) * pow(rangej, 2)) / (4.0 * M_PI * pj * pow(10, gj / 10) * pow(10, gprime / 10) * b * pow(10, loss / 10)), 0.25) / 1000.0;
}
