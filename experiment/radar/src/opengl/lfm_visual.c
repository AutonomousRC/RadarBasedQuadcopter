#include "lfm_visual.h"
#include "ogl_helper.h"

void lfm_visual(void)
{
	static char label[100];

	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glColor3f(0, 1, 0);

	glBegin(GL_LINE_LOOP);
	glVertex3f(130.0, -80.0, 0.0);
	glVertex3f(-130.0, -80.0, 0.0);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(130.0, 80.0, 0.0);
	glVertex3f(-130.0, 80.0, 0.0);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(-130.0, 80.0, 0.0);
	glVertex3f(-130.0, -80.0, 0.0);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(130.0, 80.0, 0.0);
	glVertex3f(130.0, -80.0, 0.0);
	glEnd();

#if 0
	glBegin(GL_LINE_LOOP);
	glVertex3f(0.0, 80.0, 0.0);
	glVertex3f(0.0, -80.0, 0.0);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(130.0, 0.0, 0.0);
	glVertex3f(-130.0, 0.0, 0.0);
	glEnd();
#endif

	glutSwapBuffers();
}
