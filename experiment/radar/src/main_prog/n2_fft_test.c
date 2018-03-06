#include <stdio.h>
#include "n2_fft.h"
#include "complex_math.h"

int main(void)
{
	int i;
	double sample[32] = {0};
	comp freq[256] = {0};

	for(i = 0; i < 32; i++)
		sample[i] = 1.0;

	n2_fft(sample, freq, 32, 256);

	for(i = 0; i < 256; i++)
	{
		printf("freq[%d].re = %lf\n", i, freq[i].re);
		printf("freq[%d].im = %lf\n", i, freq[i].im);
	}

	return 0;
}
