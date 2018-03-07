//#include "ogl_helper.h"
#include "init_ogl.h"
//#include "math_tech.h"
#include "sinc_function.h"
#include <stdio.h>
//#include <math.h>

int main(int argc, char **argv)
{
	char str[64] = "3D Sinc Function";

	ogl_3d_init(argc, argv, str);
	ogl_shader_process(sinc_function);

	return 0;
}
