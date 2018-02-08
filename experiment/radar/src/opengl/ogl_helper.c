#include "radar_basic.h"
#include "ogl_helper.h"
#include <string.h>

double snr[3][1000];
double rcs[4];
double range[1001];

float v_ratio, h_ratio;
float v_interval, h_interval;

// for convert to plotting
double snr_conv[3][1000];
double range_conv[1001];

void drawString(char *s)
{
	unsigned int i;

	for(i = 0; i < strlen(s); i++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, s[i]);
}

void drawStringBig(char *s)
{
	unsigned int i;

	for(i = 0; i < strlen(s); i++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, s[i]);
}

void slc(float start, float end, int num, float *arr)
{
	int i;
	float interval = (end - start) / num;

	for(i = 0; i < num; i++)
		arr[i] = start + interval * (i + 1);
}

void p_arr(double *arr, int num)
{
	int i;

	for(i = 0; i < num; i++)
		printf("arr[%d] = %lf\n", i, arr[i]);
}

void draw_rcs_dr_curve(void)
{
	int i, j;

	for(i = 0; i < 3; i++)
	{
		if(i == 0)
			glColor3f(89.0/255.0, 116.0/255.0, 238.0/255.0);
		else if(i == 1)
			glColor3f(135.0/255.0, 117.0/255.0, 204.0/255.0);
		else if(i == 2)
			glColor3f(240.0/255.0, 90.0/255.0, 53.0/255.0);

		glBegin(GL_LINES);

		for(j = 0; j < 1000; j++)
		{
			glVertex2f(range_conv[j], snr_conv[i][j]);
			//glVertex2f(range[j] - 150.0 + (j * v_interval), snr[i][j] * h_ratio - 80.0 + (j * h_interval));
			//glVertex2f(range[j] * v_ratio, snr[i][j] * h_ratio);
			//printf("range[%d] = %lf, snr[%d][%d] = %lf\n", j, range[j], i, j, snr[i][j]);
		}

		glEnd();
	}
}

float get_ratio(float gs, float ge, float os, float oe)
{
	return (ge - gs) / (oe - os);
}

float get_interval(float s, float e, float n)
{
	return (e - s) / n;
}

void convert_orig2graph(void)
{
	int i, j;
	double x_int, y_int;
	double x_ratio, y_ratio;
	double tmp_r[1001] = {0};
	double tmp_snr[3][1000] = {0};

	// range: 20 ~ 170 km
	// snr: -8.x ~ 29.x dB
	// graph x: -130 ~ 130
	// graph y: -80 ~ 80

	/* Process

	   1. Make sure zero criteria.
	      0 ~ 260:    0 - a, 260 - b
	      0 ~ 150 km: 0 - c, 150 - d

	   2. ratio = b / d = 260 / 150

	   3. interval = (d - c) / 999

	   4. index = 0 ~ 999: total 1000
	      index * interval * ratio = 0
	      index * interval * ratio = 260

	   5. 0 - 130 = -130
	      260 - 130 = 130 */

	for(i = 0; i < 1000; i++)
		tmp_r[i] = range[i] - range[0];

	x_ratio = 260 / tmp_r[i - 1];
	x_int = tmp_r[i - 1] / 999;

	for(i = 0; i < 1000; i++)
		range_conv[i] = i * x_ratio * x_int - 130.0;

	/* SNR part Process

	   -80 ~ 80 to -10 ~ 50

	   1. Make sure zero criteria.
	      0 ~ 160:    0 - a, 160 - b
	      0 ~ 60 km:  0 - c, 60 - d

	   2. ratio = b / d = 160 / 60

	   3. interval = (d - c) / 999

	   4. index = 0 ~ 999: total 1000
	      index * interval * ratio = 0
	      index * interval * ratio = 160

	   5. 0 - 80 = -80
	      160 - 80 = 80 */

	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 1000; j++)
			tmp_snr[i][j] = snr[i][j] + 10.0;

		y_ratio = 160 / 60;

		for(j = 0; j < 1000; j++)
			snr_conv[i][j] = y_ratio * tmp_snr[i][j] - 80.0;
	}
}

void display(void)
{
	static char label[100];
	float vert[16] = {0};
	float hori[7] = {0};
	int start_val = 20;
	int hori_val = -10;
	int i;

	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glColor3f(0, 1, 0);

	glBegin(GL_LINE_LOOP);
	glVertex3f(130.0, -80.0, 0.0);
	glVertex3f(-130.0, -80.0, 0.0);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(-130.0, 80.0, 0.0);
	glVertex3f(-130.0, -80.0, 0.0);
	glEnd();

	// vertical - 15, horizontal - 6
	slc(-130, 130, 15, vert);
	slc(-80, 80, 6, hori);

	v_ratio = get_ratio(-130, 130, 20, 170);
	h_ratio = get_ratio(-80, 80, -10, 50);

	//printf("v_ratio = %f, h_ratio = %f\n", v_ratio, h_ratio);

	v_interval = get_interval(-130, 130, 1000);
	h_interval = get_interval(-80, 80, 1000);

	//printf("v_interval = %f, h_interval = %f\n", v_interval, h_interval);

	convert_orig2graph();

#if 0
	p_arr(range_conv, 1000);

	for(i = 0; i < 3; i++)
		p_arr(snr_conv[i], 1000);

	p_arr(vert, 15);
	p_arr(hori, 6);
#endif

	for(i = 0; i < 15; i++)
	{
		glBegin(GL_LINE_LOOP);
		glVertex3f(vert[i], 80.0, 0.0);
		glVertex3f(vert[i], -80.0, 0.0);
		glEnd();
	}

	for(i = 0; i < 6; i++)
	{
		glBegin(GL_LINE_LOOP);
		glVertex3f(130.0, hori[i], 0.0);
		glVertex3f(-130.0, hori[i], 0.0);
		glEnd();
	}

	glColor3f(1.0, 1.0, 1.0);
	sprintf(label, "SNR");
	glRasterPos2f(-145.0, -1.5);
	drawStringBig(label);

	sprintf(label, "Detection Range");
	glRasterPos2f(-15.0, -95);
	drawStringBig(label);

	glColor3f(0.0, 1.0, 1.0);

	for(i = 0; i < 16; i++)
	{
		sprintf(label, "%d", start_val);

		start_val += 10;

		if(i)
			glRasterPos2f(vert[i - 1] - 3, -88);
		else
			glRasterPos2f(-133, -88);

		drawStringBig(label);
	}

	for(i = 0; i < 7; i++)
	{
		sprintf(label, "%d", hori_val);

		hori_val += 10;

		if(i)
			glRasterPos2f(-129.0, hori[i - 1] + 1);
		else
			glRasterPos2f(-129.0, -79.0);

		drawStringBig(label);
	}

	glColor3f(135.0/255.0, 117.0/255.0, 204.0/255.0);
	sprintf(label, "RCS = 1.0");
	glRasterPos2f(-127.2, 40);
	drawStringBig(label);

	glColor3f(89.0/255.0, 116.0/255.0, 238.0/255.0);
	sprintf(label, "RCS = 0.1");
	glRasterPos2f(-127.2, 15.5);
	drawStringBig(label);

	glColor3f(240.0/255.0, 90.0/255.0, 53.0/255.0);
	sprintf(label, "RCS = 0.01");
	glRasterPos2f(-127.2, -5.0);
	drawStringBig(label);

	draw_rcs_dr_curve();

	glutSwapBuffers();
}

void reshape(int w, int h)
{
	GLfloat n_range = 100.0f;

	if(h == 0)
		h = 1;

	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	if(w <= h)
		glOrtho(-n_range, n_range, -n_range * h / w, n_range * h / w, -n_range, n_range);
        else
		glOrtho(-n_range * w / h, n_range * w / h, -n_range, n_range, -n_range, n_range);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}
