#ifndef __COMPLEX_MATH_H__
#define __COMPLEX_MATH_H__

typedef struct
{
	double re;
	double im;
} complex;

typedef struct
{
	char re;
	char im;
} char_complex;

void euler_formula(double *, complex *, int);
void complex_abs(complex *, double *, int);
void complex_sqrt(complex *);
void complex_div(complex *, complex *, double *, int);
void print_complex(complex *, int);

#endif
