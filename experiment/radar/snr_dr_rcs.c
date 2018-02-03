#include "radar_basic.h"
#include "radar_range_equation.h"
#include "ogl_helper.h"

void slice(float start, float end, int num, float *arr)
{
	int i;
	float interval = (end - start) / (num - 1);

	arr[0] = start;

	for(i = 1; i < num; i++)
		arr[i] = start + interval * i;
}

void print_arr(float *arr)
{
	int i;

	for(i = 0; arr[i] != 0; i++)
		printf("arr[%d] = %f\n", i, arr[i]);
}

int main(void)
{
	int i;
	// for radar_range_equation()
	float snr[3][1000] = {0};
	float rcs[3] = {0.1, 1, 0.01};
	float range[1001] = {0};

	slice(25, 165, 1000, range);

	print_arr(range);

#if 0
	for(i = 0; i < 3; i++)
		radar_range_equation(1500000, 5600000000, 45, rcs[i], 290, 5000000, 3, 6, 60.78, &snr);
#endif

	return 0;
}
