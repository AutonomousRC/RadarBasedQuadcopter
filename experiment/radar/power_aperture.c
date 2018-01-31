#include "radar_basic.h"
#include "math_tech.h"

void power_aperture(float snr, float tsc, float sigma, float range, float nf,
                float to, float loss, float az_angle, float el_angle, float *pap)
{
	float tsc_db, sigma_db, four_pi_db, k_db, to_db, range_pwr4_db, omega_db;
	tsc_db = 10 * log10(tsc);
	sigma_db = 10 * log10(sigma);
	four_pi_db = 10 * log10(4.0 * M_PI);
	k_db = 10 * log10(1.38 * pow(10, -23));
	to_db = 10 * log10(to);
	range_pwr4_db = 10 * log10(pow(range, 4));
	// 57.296 = 1 rad
	omega_db = 10 * log10((az_angle / 57.296) * (el_angle / 57.296));

	*pap = snr + four_pi_db + k_db + to_db + nf + loss + range_pwr4_db + omega_db - sigma_db - tsc_db;
}
