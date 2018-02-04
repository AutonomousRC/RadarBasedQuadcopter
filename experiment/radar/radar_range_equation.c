#include "radar_basic.h"
#include "math_tech.h"

void radar_range_equation(float pt, float freq, float g, float sigma, float to,
			float b, float nf, float loss, float range, float *snr)
{
	float lambda, p_peak, lambda_sqdb, sigmadb, four_pi_cub;
	float k_db, to_db, b_db, range_pwr4_db, num, den;
	float tmp;

	lambda = get_wavelength(freq);
	p_peak = 10 * log10(pt);
	lambda_sqdb = 10 * log10(pow(lambda, 2));
	sigmadb = 10 * log10(sigma);
	four_pi_cub = 10 * log10(pow((4.0 * M_PI), 3));
	k_db = 10 * log10(1.38 / pow(10, 23));
	to_db = 10 * log10(to);
	b_db = 10 * log10(b);
	range_pwr4_db = 10 * log10(pow(range * 1000, 4));
	num = p_peak + g + g + lambda_sqdb + sigmadb;
	//den = four_pi_cub + k_db + to_db + b_db + nf + loss + range_pwr4_db;
	den = four_pi_cub + k_db + to_db + b_db + nf + loss;

	printf("lambda = %f \n", lambda);
	printf("p_peak = %f dB\n", p_peak);
	printf("lambda_sqdb = %f dB\n", lambda_sqdb);
	printf("sigmadb = %f dB\n", sigmadb);
	printf("four_pi_cub = %f dB\n", four_pi_cub);
	printf("k_db = %f dB\n", k_db);
	printf("to_db = %f dB\n", to_db);
	printf("b_db = %f dB\n", b_db);
	printf("range_pwr4_db = %f dB\n", range_pwr4_db);

	tmp = num - den - 20;
	printf("kBT_o = %f\n", k_db + to_db + b_db);
	printf("tmp = %f\n", tmp);
	printf("range^4 = %f\n", pow(10, tmp/10));
	printf("range = %f\n", pow(pow(10, tmp/10), 0.25) / 1000);

	den = four_pi_cub + k_db + to_db + b_db + nf + loss + range_pwr4_db;

	printf("snr = %f\n", p_peak + 2 * g + lambda_sqdb + sigmadb - four_pi_cub - k_db - to_db - b_db - nf - loss - range_pwr4_db);

	/* 2018/01/30 - It makes some confusion.
	   X dB = 10 * log_10(R^4 * 1000)
	   (X / 40) dB = log_10(R * 1000) */

	*snr = num - den;
}

void rre_vec(double pt, double freq, double g, double sigma, double to,
	double b, double nf, double loss, double *range, double *snr)
{
	int i;
	double tmp;
	double lambda, p_peak, lambda_sqdb, sigmadb, four_pi_cub;
	double k_db, to_db, b_db, num, den;
	double range_pwr4_db[1000] = {0};

	lambda = get_wavelength(freq);
	p_peak = 10 * log10(pt);
	lambda_sqdb = 10 * log10(pow(lambda, 2));
	sigmadb = 10 * log10(sigma);
	four_pi_cub = 10 * log10(pow((4.0 * M_PI), 3));
	k_db = 10 * log10(1.38 / pow(10, 23));
	to_db = 10 * log10(to);
	b_db = 10 * log10(b);

	for(i = 0; i < 1000; i++)
		range_pwr4_db[i] = 10 * log10(pow(range[i] * 1000, 4));

	num = p_peak + g + g + lambda_sqdb + sigmadb;
	den = four_pi_cub + k_db + to_db + b_db + nf + loss;

	for(i = 0; i < 1000; i++)
	{
		tmp = den;
		tmp += range_pwr4_db[i];
		snr[i] = num - tmp;
	}
}
