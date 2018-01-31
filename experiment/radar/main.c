#include "radar_basic.h"
#include "pulse_train.h"
#include "range_resolution.h"
#include "doppler_frequency.h"
#include "radar_range_equation.h"
#include "power_aperture.h"
#include "self_screening_jammer.h"

int main(void)
{
	int ret;
	// for pulse_train()
	float dt, prf, pav, ep, ru;
	// for range_resolution()
	float delta_R;
	// for doppler_frequency()
	float df, td;
	// for radar_range_equation()
	float snr;
	// for power_aperture()
	float pap;
	// for self_screening_jammer()
	float br_range;

	// test for airborne pulsed radar(peak power 10 KW)
	pulse_train(0.000015, 0.0001, 10000, &dt, &prf, &pav, &ep, &ru);

	printf("dt = %f\n", dt);
	printf("prf = %f\n", prf);
	printf("pav = %f\n", pav);
	printf("ep = %f\n", ep);
	printf("ru = %f\n", ru);

	// unambiguous range of 100 km, bandwidth 0.5 MHz
	if(ret = range_resolution(500000, &delta_R, 1))
	{
		printf("range_resolution() error!\n");
		exit(-1);
	}

	printf("delta_R = %f\n", delta_R);

	// test for X-Band 10 GHz Radar with 10 Mach number
	if(ret = doppler_frequency(10000000000.0, 0.0, 250.0, -175.0, &df, &td))
	{
		printf("doppler_frequency() error!\n");
		exit(-1);
	}

	printf("df = %f\n", df);
	printf("td = %f\n", td);

	// radar range equation test
	// Peak Power = 1.5 MW, Operating Frequency = 5.6 GHz, Antenna Gain = 45 dB,
	// Effective Temperature = 290 K, Pulse Width = 0.2 us.
	// We can detect snr = 20 with max range 60.78 km.
	radar_range_equation(1500000, 5600000000, 45, 0.1, 290, 5000000, 3, 6, 60.78, &snr);

	printf("snr = %f dB\n", snr);

	// Scan Time = 2 s, Noise Figure = 8 dB, Losses = 6 dB, Search Volume = 7.4 Steradians
	// Range of Interest = 75 km, Required SNR = 20 dB, Target Cross Section = 3.162 m^2
	// Effective Temperature = 290 K, Search Volume Azimuth = 180 degree,
	// Search Volume Elevation = 135 degree.
	power_aperture(20, 2, 3.162, 75000, 8, 290, 6, 180, 135, &pap);

	printf("pap = %f dB\n", pap);

	// Radar Peak Power = 50 KW, Jammer Peak Power = 200 W,
	// Radar Operating Bandwidth = 667 KHz, Jammer Bandwidth = 50 MHz,
	// Radar & Jammer Losses = 0.1 dB, Target Cross Section = 10 m^2,
	// Radar Antenna Gain = 35 dB, Jammer Antenna Gain = 10 dB,
	// Radar Operating Frequency = 5.6 GHz
	self_screening_jammer(50000, 35, 5600000000, 10, 667000, 0.1, 200, 50000000, 10, 0.1, &br_range);

	printf("br_range(cross-over range) = %f km\n", br_range);

	return 0;
}
