#ifndef __PULSE_TRAIN_H__
#define __PULSE_TRAIN_H__

#define RT_PI		3.14159265358979323846

float angle2radian(float);
float get_wavelength(float);
void slice_section(double, double, double, double *);
void linear_slice(double, double, int, double *);
void pow_vec(double *, double *, int, double);
void p_arr(double *, int);
void swap(double *, int, int);
int partition(double *, int, int, int);
void quicksort(double *, int, int);
double find_max(double *, int);
double find_min(double *, int);
void divide_vec(double *, double, int);
void volt_base_log(double *, double, int);
double dist_2d(double, double);

#endif
