#include "radar_basic.h"
#include "radar_range_equation.h"
#include "ogl_helper.h"

// for radar_range_equation()
extern double snr[3][1000];
extern double rcs[4];
extern double range[1001];

void slice(double start, double end, int num, double *arr)
{
	int i;
	double interval = (end - start) / (num - 1);

	arr[0] = start;

	for(i = 1; i < num; i++)
		arr[i] = start + interval * i;
}

void print_arr(double *arr)
{
	int i;

	for(i = 0; arr[i] != 0; i++)
		printf("arr[%d] = %lf\n", i, arr[i]);
}

void print_darr(double (*arr)[1000])
{
	int i, j;

	for(i = 0; i < 3; i++)
		for(j = 0; j < 1000; j++)
			printf("arr[%d][%d] = %.24lf\n", i, j, arr[i][j]);
}

int main(int argc, char **argv)
{
	int i;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE);
	glutInitWindowSize(1200, 800);
	//glutInitWindowPosition(0, 0);
	glutCreateWindow("SNR vs Detection Range");

	slice(20, 180, 1000, range);

	rcs[0] = 0.1;
	rcs[1] = 1.0;
	rcs[2] = 0.01;

#if DEBUG
	print_arr(range);
#endif

	// radar range equation
	for(i = 0; i < 3; i++)
		rre_vec(1500000, 5600000000, 45, rcs[i], 290, 5000000, 3, 6, range, snr[i]);

#if DEBUG
	print_darr(snr);
#endif

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMainLoop();

	return 0;
}
