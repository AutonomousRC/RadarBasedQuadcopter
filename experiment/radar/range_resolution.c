#include "radar_basic.h"

int range_resolution(float variable, float *delta_R, int sel)
{
	int nread;
	int in_num;
	char buf[32] = {0};

#if 0
	printf("Enter 1 for var == Bandwidth, 2 for var == Pulse Width:\n");
	nread = read(0, buf, sizeof(buf));

	in_num = atoi(buf);

	printf("in_num = %d\n", in_num);

	if(in_num != 1 && in_num != 2)
	{
		printf("You entered wrong number\n");
		return -1;
	}
#endif

	switch(sel)
	{
		case 1:
			*delta_R = LIGHT_SPEED / (2.0 * variable);
			break;
		case 2:
			*delta_R = (LIGHT_SPEED * variable) / 2.0;
			break;
	}

	return 0;
}
