#include <stdio.h>
#include "init_require_data.h"

int main(void)
{
	double signal[800] = {0};
	init_require_data(signal, 800);

	return 0;
}
