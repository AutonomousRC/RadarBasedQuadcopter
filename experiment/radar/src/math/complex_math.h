#ifndef __COMPLEX_MATH_H__
#define __COMPLEX_MATH_H__

typedef struct
{
	double re;
	double im;
} complex;

void euler_formula(double *, complex *, int);
void complex_abs(complex *, double *, int);
void print_complex(complex *, int);

#endif
