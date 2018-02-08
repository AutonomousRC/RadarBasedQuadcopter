#include "radar_basic.h"
#include "math_tech.h"

void self_screening_jammer(float pt, float g, float freq, float sigma, float br,
                        float loss, float pj, float bj, float gj, float lossj, float *br_range)
{
	int i;
	float rv, r_db;
	float s[1000] = {0};
	float ssj[1000] = {0};
	float lambda_sqdb, sigma_db, pt_db, b_db, bj_db, pj_db, factor, s_at_br;

	lambda_sqdb = 10 * log10(pow(get_wavelength(freq), 2));

	if(loss == 0.0)
		loss = 0.000001;

	if(lossj == 0.0)
		lossj = 0.000001;

	sigma_db = 10 * log10(sigma);
	pt_db = 10 * log10(pt);
	b_db = 10 * log10(br);
	bj_db = 10 * log10(bj);
	pj_db = 10 * log10(pj);
	factor = 10 * log10(4.0 * M_PI);

	*br_range = sqrt(((pt * pow(10, g / 10)) * sigma * bj * pow(10, lossj / 10)) / (4.0 * M_PI * pj * pow(10, gj / 10) * br * pow(10, loss / 10))) / 1000.0;

	/* TODO: Add OpenGL to build graph */

	s_at_br = pt_db + 2.0 * g + lambda_sqdb + sigma_db - 3.0 * factor - 4 * 10 * log10(*br_range) - loss;

	for(i = 0, rv = 0.1; i < 1000; i++, rv += 10)
	{
		r_db = 10 * log10(rv * 1000.0);
		ssj[i] = pj_db + gj + lambda_sqdb + g + b_db - 2.0 * factor - 2.0 * r_db - bj_db - lossj + s_at_br;
		s[i] = pt_db + 2.0 * g + lambda_sqdb + sigma_db - 3.0 * factor - 4.0 * r_db - loss + s_at_br;
	}

	
}
