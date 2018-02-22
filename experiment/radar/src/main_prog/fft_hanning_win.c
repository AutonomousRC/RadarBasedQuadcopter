#include "ogl_helper.h"
#include "init_ogl.h"
#include "hanning_window.h"
#include "math_tech.h"
#include <stdio.h>
#include <math.h>

int main(int argc, char **argv)
{
	char str[64] = "Hanning Window Function";

	ogl_init(argc, argv, str);
	ogl_post_process(hanning_window);

	return 0;
}
