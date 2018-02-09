#ifndef __PULSE_TRAIN_H__
#define __PULSE_TRAIN_H__

float angle2radian(float);
float get_wavelength(float);
void slice_section(double, double, double, double *);
void linear_slice(double, double, int, double *);
void pow_vec(double *, double *, int, double);

#endif
