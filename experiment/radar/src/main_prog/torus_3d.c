#include "ogl_helper.h"
#include "init_ogl.h"
#include "draw_torus.h"

int main(int argc, char **argv)
{
	char str[25] = "Draw 3D Torus";

	ogl_3d_init(argc, argv, str);
	ogl_post_process(draw_torus);

	return 0;
}
