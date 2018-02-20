#include "ogl_helper.h"
#include "init_ogl.h"

void ogl_init(int argc, char **argv, char *str)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE);
	glutInitWindowSize(1200, 800);
	glutCreateWindow(str);
}

void ogl_3d_init(int argc, char **argv, char *str)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(1200, 800);
	glutCreateWindow(str);
}

void setup_rc(void)
{
	glClearColor(0.0f, 0.0f, 0.5f, 1.0f);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
}

void ogl_post_process(void (*disp_func)(void))
{
	glutDisplayFunc(disp_func);
	glutReshapeFunc(reshape);
	glutMainLoop();
}

void ogl_3d_post_process(void (*disp_func)(void), void (*reshape_func)(int, int), void (*timer_func)(int))
{
	glutDisplayFunc(disp_func);
	glutReshapeFunc(reshape_func);
	setup_rc();

	if(timer_func)
		glutTimerFunc(33, timer_func, 1);

	glutMainLoop();
}
