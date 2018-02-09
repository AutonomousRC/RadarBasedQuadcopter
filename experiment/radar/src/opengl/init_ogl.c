#include "ogl_helper.h"
#include "init_ogl.h"

void ogl_init(int argc, char **argv, char *str)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE);
	glutInitWindowSize(1200, 800);
	glutCreateWindow(str);
}

void ogl_post_process(void (*disp_func)(void))
{
	glutDisplayFunc(disp_func);
	glutReshapeFunc(reshape);
	glutMainLoop();
}
