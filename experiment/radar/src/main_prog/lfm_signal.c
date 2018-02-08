#include "ogl_helper.h"
#include "lfm_visual.h"
#include "math_tech.h"
#include <stdio.h>
#include <math.h>

int main(int argc, char **argv)
{
	double lfm_bandwidth = 15 * pow(10, 6);
	double tau = pow(10, -6);
	double ts = pow(10, -9);
	double beta = lfm_bandwidth / tau;

	double t[1002] = {0};

	// t = 0 ~ 0.000001 each slice 0.000000001
	slice_section(0.0, tau, ts, t);
	//printf("check = %lf\n", tau / ts);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE);
	glutInitWindowSize(1200, 800);
	glutCreateWindow("Typical LFM Waveforms");

	glutDisplayFunc(lfm_visual);
	glutReshapeFunc(reshape);
	glutMainLoop();

	return 0;
}
