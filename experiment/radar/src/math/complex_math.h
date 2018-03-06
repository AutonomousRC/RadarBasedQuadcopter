#ifndef __COMPLEX_MATH_H__
#define __COMPLEX_MATH_H__

typedef struct
{
	double re;
	double im;
} comp;

typedef struct
{
	char re;
	char im;
} char_complex;

void euler_formula(double *, comp *, int);
void complex_abs(comp *, double *, int);
void complex_sqrt(comp *);
void complex_div(comp *, comp *, double *, int);
void print_complex(comp *, int);

#endif
