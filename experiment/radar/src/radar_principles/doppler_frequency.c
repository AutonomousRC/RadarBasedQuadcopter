#include "radar_basic.h"
#include "math_tech.h"

int doppler_frequency(float freq, float ang, float radar_vel, float tar_vel,
			float *doppler_freq, float *time_dil)
{
	int nread;
	int in_num;
	char buf[32] = {0};
	float ang_rad, lambda, vel;

	printf("Enter 1 for closing target, 2 for opening target:\n");
	nread = read(0, buf, sizeof(buf));

	in_num = atoi(buf);

	printf("in_num = %d\n", in_num);

	if(in_num != 1 && in_num != 2)
	{
		printf("You entered wrong number\n");
		return -1;
	}
	
	ang_rad = angle2radian(ang);
	lambda = get_wavelength(freq);

	printf("ang_rad = %f\n", ang_rad);
	printf("lambda = %f\n", lambda);
	printf("cos(ang_rad) = %f\n", cos(ang_rad));

	switch(in_num)
	{
		case 1:
			vel = radar_vel - tar_vel;
			*doppler_freq = 2.0 * vel * cos(ang_rad) / lambda;
			*time_dil = (LIGHT_SPEED - vel) / (LIGHT_SPEED + vel);
			break;
		case 2:
			vel = radar_vel + tar_vel;
			*doppler_freq = 2.0 * vel * cos(ang_rad) / lambda;
			*time_dil = (LIGHT_SPEED + vel) / (LIGHT_SPEED - vel);
			break;
	}

	return 0;
}
