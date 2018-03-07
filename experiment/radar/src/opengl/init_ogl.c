#include <stdio.h>
#include <stdlib.h>
#include <GL/glew.h>

#include "ogl_helper.h"
#include "init_ogl.h"
#include "sinc_function.h"

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

void ogl_shader_process(void (*disp_func)(void))
{
    GLenum glew_status = glewInit();

    if(GLEW_OK != glew_status)
    {
        fprintf(stderr, "error: %s\n", glewGetErrorString(glew_status));
        exit(1);
    }

    if(!GLEW_VERSION_2_0)
    {
        fprintf(stderr, "No support for OpenGL 2.0 found\n");
        exit(1);
    }

    GLint max_units;

    glGetIntegerv(GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS, &max_units);

    if(max_units < 1)
    {
        fprintf(stderr, "GPU doesn't have any vertex texture image units\n");
        exit(1);
    }

    if(init_shaders())
    {
        glutDisplayFunc(disp_func);
        glutIdleFunc(disp_func);
        glutMainLoop();
    }

    free_resources();
}
