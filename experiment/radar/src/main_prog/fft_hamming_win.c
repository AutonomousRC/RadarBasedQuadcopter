#include "ogl_helper.h"
#include "init_ogl.h"
#include "hamming_window.h"
#include "math_tech.h"
#include <stdio.h>
#include <math.h>

int main(int argc, char **argv)
{
	char str[64] = "Hamming Window Function";

	ogl_init(argc, argv, str);
	ogl_post_process(hamming_window);

	return 0;
}
